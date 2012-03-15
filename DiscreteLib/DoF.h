
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
 * \brief Mapping of DoFs for one variable
 *
 *  - discrete point id <-> equation id
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
        return (_list_ghost_node_id.count(discrete_pt_id) > 0);
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

    void getListOfGhostNodeID(std::vector<size_t> &list_ghost)
    {
        list_ghost.assign(_list_ghost_node_id.begin(), _list_ghost_node_id.end());
    }

    size_t getOrder() const {return _order;};
private:
    std::vector<long> _map_node_id2eqs_id;
    size_t _order;
    std::set<size_t> _list_inactive_node_id;
    std::set<size_t> _list_ghost_node_id;

    DISALLOW_COPY_AND_ASSIGN(DofMap);
};


}
