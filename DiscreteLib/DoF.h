
#pragma once

#include <vector>
#include <map>
#include <set>
#ifdef USE_MPI
#include "mpi.h"
#endif

#include "Base/BidirectionalMap.h"
#include "Base/CodingTools.h"
#include "DiscreteVector.h"

namespace DiscreteLib
{

/**
 * \brief Mapping of DOF id and mesh node id
 *
 *  - discrete point id <-> eqs id
 */
class DofMap
{
public:
    DofMap(size_t discrete_points_size, size_t order=1) : _map_node_id2eqs_id(discrete_points_size), _order(order)
    {
    }
    DofMap(size_t discrete_points_size, std::set<size_t>* list_inactive_node_id, size_t order=1) : _map_node_id2eqs_id(discrete_points_size), _order(order)
    {
        if (list_inactive_node_id!=0)
            _list_inactive_node_id.insert(list_inactive_node_id->begin(), list_inactive_node_id->end());
    }
    DofMap(size_t discrete_points_size, std::set<size_t>* list_ghost_node_id, std::set<size_t>* list_inactive_node_id, size_t order=1) : _map_node_id2eqs_id(discrete_points_size), _order(order)
    {
        if (list_inactive_node_id!=0)
            _list_inactive_node_id.insert(list_inactive_node_id->begin(), list_inactive_node_id->end());
        if (list_ghost_node_id!=0)
            _list_ghost_node_id.insert(list_ghost_node_id->begin(), list_ghost_node_id->end());
    }

    void getListOfEqsID( const std::vector<size_t>& vec_pt_id, std::vector<long>& vec_eqs_id ) const
    {
        getListOfEqsID(vec_pt_id, vec_pt_id.size(), vec_eqs_id);
    }

    void getListOfEqsID( const std::vector<size_t>& vec_pt_id, size_t n, std::vector<long>& vec_eqs_id ) const
    {
        vec_eqs_id.resize(n);
        for (size_t i=0; i<n; i++) {
            vec_eqs_id[i] = getEqsID(vec_pt_id[i]);
        }
    }

    bool isActiveDoF(size_t discrete_pt_id) const 
    {
        return (_list_inactive_node_id.count(discrete_pt_id) == 0);
    }

    bool isGhostDoF(size_t discrete_pt_id) const 
    {
        return (_list_ghost_node_id.count(discrete_pt_id) == 0);
    }

    long getEqsID(size_t discrete_pt_id) const
    {
        return _map_node_id2eqs_id[discrete_pt_id];
    }

    void setEqsID(size_t discrete_pt_id, long eqs_id)
    {
        _map_node_id2eqs_id[discrete_pt_id] = eqs_id;
    }

    size_t setEqsIDSequnetual(size_t eqs_id_begin);

    size_t getNumberOfDiscretePoints() const 
    {
        return _map_node_id2eqs_id.size();
    }

    size_t getNumberOfGhostPoints() const 
    {
        return _list_ghost_node_id.size();
    };

    size_t getNumberOfActiveDoFs() const 
    {
        return _map_node_id2eqs_id.size() - _list_inactive_node_id.size();
    }

    size_t getNumberOfActiveDoFsWithoutGhost() const
    {
        size_t count = 0;
        for (size_t i=0; i<_map_node_id2eqs_id.size(); i++) {
            if (isActiveDoF(i) && !isGhostDoF(i) ) count++;
        }
        return count;
    }

    size_t getOrder() const {return _order;};
private:
    std::vector<long> _map_node_id2eqs_id;
    size_t _order;
    std::set<size_t> _list_inactive_node_id;
    std::set<size_t> _list_ghost_node_id;

    DISALLOW_COPY_AND_ASSIGN(DofMap);
};

/**
 * \brief Dof map manger
 */
class DofMapManager
{
public:
    enum NumberingType
    {
        BY_DOF,
        BY_POINT
    };

    DofMapManager() : _total_pt(0), _total_dofs(0) {};
    virtual ~DofMapManager()
    {
        Base::releaseObjectsInStdVector(_map_var2dof);
    }

    size_t addDoF(size_t discrete_points_size, size_t order=1, size_t mesh_id=0);

    size_t addDoF(size_t discrete_points_size, std::set<size_t>* list_inactive_node_id, size_t order=1, size_t mesh_id=0);

    size_t addDoF(size_t discrete_points_size, std::set<size_t>* ghost_nodes, std::set<size_t>* list_inactive_node_id, size_t order=1, size_t mesh_id=0);

    virtual void construct(NumberingType num=BY_DOF, size_t offset=0);

    size_t getNumberOfDof() const
    {
        return _map_var2dof.size();
    }

    const DofMap* getDofMap(size_t var_id) const
    {
        return _map_var2dof.at(var_id);
    }

    DofMap* getDofMap(size_t var_id)
    {
        return _map_var2dof.at(var_id);
    }

    size_t getTotalNumberOfDiscretePoints() const
    {
        return _total_pt;
    }

    size_t getTotalNumberOfActiveDoFs() const
    {
        return _total_dofs;
    }

    size_t getTotalNumberOfActiveDoFsWithoutGhost() const
    {
        size_t count = 0;
        for (size_t i=0; i<_map_var2dof.size(); i++) {
            count += _map_var2dof[i]->getNumberOfActiveDoFsWithoutGhost();
        }
        return count;
    }


    void getListOfEqsID(const std::vector<size_t> &ele_node_ids, const std::vector<size_t> &ele_node_size_order, std::vector<long> &local_dofmap) const;

    size_t getTotalNumberOfDOFs(const std::vector<size_t> &list_vec_entry_id, const std::vector<size_t> &list_vec_size_for_order) const;

    /// create a subset of vector u corresponding to the given vector index
    void getLocalVector(const std::vector<size_t> &list_vec_entry_id, const std::vector<size_t> &list_vec_size_for_order, const std::vector<DiscreteVector<double>*> &list_multiple_u, std::vector<double> &local_u) const;

protected:
    void setTotalDoFs(size_t n) {_total_dofs = n;};
    std::map<size_t, std::vector<DofMap*> >& getMapMsh2Dof() {return _map_msh2dof;}
private:
    std::vector<DofMap*> _map_var2dof;
    std::map<size_t, std::vector<DofMap*> > _map_msh2dof;
    size_t _total_pt;
    size_t _total_dofs;

    DISALLOW_COPY_AND_ASSIGN(DofMapManager);
};

#ifdef USE_MPI

class MPIDofMapManager : public DofMapManager
{
public:
    MPIDofMapManager(MPI_Comm comm, Base::BidirectionalMap<size_t, size_t> &map_global_local) : _map_global_local(map_global_local) 
    {
        _my_comm = comm;
        MPI_Comm_size(_my_comm,&_nprocs);
        MPI_Comm_rank(_my_comm,&_my_rank);
    };

    virtual ~MPIDofMapManager()
    {
    }

    void construct()
    {
        size_t dof_end_prev_rank = 0;
        //receive
        MPI_Status status;
        if (_my_rank>0)
            MPI_Recv(&dof_end_prev_rank, 1, MPI_INT, _my_rank-1, MPI_ANY_TAG, _my_comm, &status);
        //set
        DofMapManager::construct(BY_POINT, dof_end_prev_rank);
        size_t dof_end_current = getTotalNumberOfActiveDoFs();
        //send
        if (_my_rank<_nprocs-1)
            MPI_Send(&dof_end_current, 1, MPI_INT, _my_rank+1, MPI_ANY_TAG, _my_comm);
    }

private:
    MPI_Comm _my_comm;
    int _my_rank;
    int _nprocs;
    Base::BidirectionalMap<size_t, size_t> _map_global_local;

    DISALLOW_COPY_AND_ASSIGN(MPIDofMapManager);
};
#endif

}
