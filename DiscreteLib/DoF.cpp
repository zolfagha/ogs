
#include "DoF.h"


namespace DiscreteLib
{

size_t DofMap::setEqsIDSequnetual(size_t eqs_id_begin)
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

size_t DofMapManager::addDoF(size_t discrete_points_size, size_t order, size_t mesh_id)
{
    return addDoF(discrete_points_size, 0, order, mesh_id);
}

size_t DofMapManager::addDoF(size_t discrete_points_size, std::set<size_t>* list_inactive_node_id, size_t order, size_t mesh_id)
{
    return addDoF(discrete_points_size, 0, list_inactive_node_id, order, mesh_id);
}

size_t DofMapManager::addDoF(size_t discrete_points_size, std::set<size_t>* ghost_nodes, std::set<size_t>* list_inactive_node_id, size_t order, size_t mesh_id)
{
    _total_pt += discrete_points_size;
    DofMap *dof = new DofMap(discrete_points_size, ghost_nodes, list_inactive_node_id, order);
    _map_var2dof.push_back(dof);
    _map_msh2dof[mesh_id].push_back(dof);
    return _map_var2dof.size()-1;
}

void DofMapManager::construct(NumberingType num, size_t offset)
{
    if (num==BY_DOF) {
        //order by dof
        size_t eqs_id = offset;
        for (size_t i=0; i<_map_var2dof.size(); i++) {
            DofMap *dof = _map_var2dof[i];
            eqs_id = dof->setEqsIDSequnetual(eqs_id);
        }
        _total_dofs = eqs_id;
    } else {
        //order by discrete points
        size_t eqs_id = offset;
        for (std::map<size_t, std::vector<DofMap*> >::iterator itr=_map_msh2dof.begin(); itr!=_map_msh2dof.end(); itr++) {
            std::vector<DofMap*> *vec = &itr->second;
            size_t pt_size = vec->at(0)->getNumberOfDiscretePoints();
            // first, real nodes
            for (size_t i=0; i<pt_size; i++) {
                for (size_t j=0; j<vec->size(); j++) {
                    DofMap *dof = vec->at(j);
                    if (dof->isGhostDoF(i)) continue;
                    if (dof->isActiveDoF(i))
                        dof->setEqsID(i, eqs_id++);
                    else
                        dof->setEqsID(i, -1);
                }
            }
            // next ghost nodes
            for (size_t i=0; i<pt_size; i++) {
                for (size_t j=0; j<vec->size(); j++) {
                    DofMap *dof = vec->at(j);
                    if (!dof->isGhostDoF(i)) continue;
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

void DofMapManager::getListOfEqsID(const std::vector<size_t> &ele_node_ids, const std::vector<size_t> &ele_node_size_order, std::vector<long> &local_dofmap) const
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

size_t DofMapManager::getTotalNumberOfDOFs(const std::vector<size_t> &list_vec_entry_id, const std::vector<size_t> &list_vec_size_for_order) const
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
void DofMapManager::getLocalVector(const std::vector<size_t> &list_vec_entry_id, const std::vector<size_t> &list_vec_size_for_order, const std::vector<DiscreteVector<double>*> &list_multiple_u, std::vector<double> &local_u) const
{
    local_u.clear();
    const size_t n_dof = getNumberOfDof();
    for (size_t i=0; i<n_dof; i++) {
        const DofMap *dofMap = getDofMap(i);
        const size_t vec_size = list_vec_size_for_order[dofMap->getOrder()-1];
        const DiscreteVector<double> &dof_u = *list_multiple_u[i];
        for (size_t j=0; j<vec_size; j++) {
            local_u.push_back(dof_u[list_vec_entry_id[j]]);
        }
    }
}

} //end
