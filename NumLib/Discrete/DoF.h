
#pragma once

#include <vector>
#include <map>
#include "Base/CodingTools.h"
#include "MeshLib/Core/IMesh.h"
#include "MeshLib/Core/IElement.h"

namespace NumLib
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
    DofMap(size_t discrete_points_size, std::set<size_t> &list_inactive_node_id, size_t order=1) : _map_node_id2eqs_id(discrete_points_size), _order(order)
    {
        _list_inactive_node_id.insert(list_inactive_node_id.begin(), list_inactive_node_id.end());
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

    long getEqsID(size_t discrete_pt_id) const
    {
        return _map_node_id2eqs_id[discrete_pt_id];
    }

    void setEqsID(size_t discrete_pt_id, long eqs_id)
    {
        _map_node_id2eqs_id[discrete_pt_id] = eqs_id;
    }

    size_t setEqsIDSequnetual(size_t eqs_id_begin)
    {
        size_t eqs_id = eqs_id_begin;
        for (size_t i=0; i<_map_node_id2eqs_id.size(); i++) {
            if (isActiveDoF(i))
                _map_node_id2eqs_id[i] = eqs_id++;
            else
                _map_node_id2eqs_id[i] = -1;
        }
        return eqs_id;
    }

    size_t getNumberOfDiscretePoints() const 
    {
        return _map_node_id2eqs_id.size();
    }

    size_t getNumberOfActiveDoFs() const 
    {
        return _map_node_id2eqs_id.size() - _list_inactive_node_id.size();
    }

    size_t getOrder() const {return _order;};
private:
    std::vector<long> _map_node_id2eqs_id;
    size_t _order;
    std::set<size_t> _list_inactive_node_id;

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

    size_t addDoF(size_t discrete_points_size, size_t order=1, size_t mesh_id=0)
    {
        _total_pt += discrete_points_size;
        DofMap *dof = new DofMap(discrete_points_size, order);
        _map_var2dof.push_back(dof);
        _map_msh2dof[mesh_id].push_back(dof);
        return _map_var2dof.size()-1;
    }

    size_t addDoF(size_t discrete_points_size, std::set<size_t> &list_inactive_node_id, size_t order=1, size_t mesh_id=0)
    {
        _total_pt += discrete_points_size;
        DofMap *dof = new DofMap(discrete_points_size, list_inactive_node_id, order);
        _map_var2dof.push_back(dof);
        _map_msh2dof[mesh_id].push_back(dof);
        return _map_var2dof.size()-1;
    }

    void construct(NumberingType num=BY_DOF)
    {
        if (num==BY_DOF) {
            //order by dof
            size_t eqs_id = 0;
            for (size_t i=0; i<_map_var2dof.size(); i++) {
                DofMap *dof = _map_var2dof[i];
                eqs_id = dof->setEqsIDSequnetual(eqs_id);
            }
            _total_dofs = eqs_id;
        } else {
            //order by discrete points
            size_t eqs_id = 0;
            for (std::map<size_t, std::vector<DofMap*> >::iterator itr=_map_msh2dof.begin(); itr!=_map_msh2dof.end(); itr++) {
                std::vector<DofMap*> *vec = &itr->second;
                size_t pt_size = vec->at(0)->getNumberOfDiscretePoints();
                for (size_t i=0; i<pt_size; i++) {
                    for (size_t j=0; j<vec->size(); j++) {
                        DofMap *dof = vec->at(j);
                        if (dof->isActiveDoF(i))
                            dof->setEqsID(i, eqs_id++);
                        else
                            dof->setEqsID(i, -1);
                    }
                }
            }
            _total_dofs = eqs_id;
        }
    }

    size_t getNumberOfDof() const
    {
        return _map_var2dof.size();
    }

    const DofMap* getDofMap(size_t var_id) const
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

    void getListOfEqsID(const std::vector<size_t> &ele_node_ids, const std::vector<size_t> &ele_node_size_order, std::vector<long> &local_dofmap) const
    {
        const size_t n_dof = getNumberOfDof();
        if (n_dof<2) {
            const DofMap *dofMap = getDofMap(0);
            dofMap->getListOfEqsID(ele_node_ids, ele_node_size_order[dofMap->getOrder()-1], local_dofmap);
        } else {
            for (size_t j=0; j<n_dof; j++) {
                std::vector<long> tmp_map;
                const DofMap *dofMap = getDofMap(j);
                dofMap->getListOfEqsID(ele_node_ids, ele_node_size_order[dofMap->getOrder()-1], tmp_map);
                local_dofmap.insert(local_dofmap.end(), tmp_map.begin(), tmp_map.end());
            }
        }
    }

    size_t getTotalNumberOfDOFs(const std::vector<size_t> &list_vec_entry_id, const std::vector<size_t> &list_vec_size_for_order) const
    {
        size_t s = 0;
        const size_t n_dof = getNumberOfDof();
        for (size_t i=0; i<n_dof; i++) {
            const DofMap *dofMap = getDofMap(i);
            s += list_vec_size_for_order[dofMap->getOrder()-1];
        }
        return s;
    }

    /// create a subset of vector u corresponding to the given vector index
    void getLocalVector(const std::vector<size_t> &list_vec_entry_id, const std::vector<size_t> &list_vec_size_for_order, const std::vector<std::vector<double>*> &list_multiple_u, std::vector<double> &local_u) const
    {
        local_u.clear();
        const size_t n_dof = getNumberOfDof();
        for (size_t i=0; i<n_dof; i++) {
            const DofMap *dofMap = getDofMap(i);
            const size_t vec_size = list_vec_size_for_order[dofMap->getOrder()-1];
            const std::vector<double> &dof_u = *list_multiple_u[i];
            for (size_t j=0; j<vec_size; j++) {
                local_u.push_back(dof_u[list_vec_entry_id[j]]);
            }
        }
    }

private:
    std::vector<DofMap*> _map_var2dof;
    std::map<size_t, std::vector<DofMap*> > _map_msh2dof;
    size_t _total_pt;
    size_t _total_dofs;

    DISALLOW_COPY_AND_ASSIGN(class DofMapManager);
};

}
