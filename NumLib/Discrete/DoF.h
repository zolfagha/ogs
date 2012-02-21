
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
private:
    std::vector<size_t> _map_node_id2eqs_id;
    size_t _order;

public:
    DofMap(size_t discrete_points_size, size_t order=1) : _map_node_id2eqs_id(discrete_points_size), _order(order)
    {
    }

    void getListOfEqsID( const std::vector<size_t>& vec_pt_id, std::vector<size_t>& vec_eqs_id ) const
    {
        getListOfEqsID(vec_pt_id, vec_pt_id.size(), vec_eqs_id);
    }

    void getListOfEqsID( const std::vector<size_t>& vec_pt_id, size_t n, std::vector<size_t>& vec_eqs_id ) const
    {
        vec_eqs_id.resize(n);
        for (size_t i=0; i<n; i++) {
            vec_eqs_id[i] = getEqsID(vec_pt_id[i]);
        }
    }

    size_t getEqsID(size_t discrete_pt_id) const
    {
        return _map_node_id2eqs_id[discrete_pt_id];
    }

    void setEqsID(size_t discrete_pt_id, size_t eqs_id)
    {
        _map_node_id2eqs_id[discrete_pt_id] = eqs_id;
    }

    void setEqsIDSequnetual(size_t eqs_id_begin)
    {
        for (size_t i=0; i<_map_node_id2eqs_id.size(); i++)
            _map_node_id2eqs_id[i] = eqs_id_begin + i;
    }

    size_t getNumberOfDiscretePoints() const 
    {
        return _map_node_id2eqs_id.size();
    }

    size_t getOrder() const {return _order;};
};

class DofMapManager
{
private:
    std::vector<DofMap*> _map_var2dof;
    std::map<size_t, std::vector<DofMap*> > _map_msh2dof;
    size_t _total_pt;

public:
    enum NumberingType
    {
        BY_DOF,
        BY_POINT
    };

    DofMapManager() {};
    virtual ~DofMapManager()
    {
        Base::releaseObjectsInStdVector(_map_var2dof);
    }

    size_t addDoF(size_t discrete_points_size, size_t order=1, size_t mesh_id=0)
    {
        DofMap *dof = new DofMap(discrete_points_size, order);
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
                dof->setEqsIDSequnetual(eqs_id);
                eqs_id += dof->getNumberOfDiscretePoints();
            }
            _total_pt = eqs_id;
        } else {
            //order by discrete points
            size_t eqs_id = 0;
            for (std::map<size_t, std::vector<DofMap*> >::iterator itr=_map_msh2dof.begin(); itr!=_map_msh2dof.end(); itr++) {
                std::vector<DofMap*> *vec = &itr->second;
                size_t pt_size = vec->at(0)->getNumberOfDiscretePoints();

                for (size_t i=0; i<pt_size; i++) {
                    for (size_t j=0; j<vec->size(); j++) {
                        DofMap *dof = vec->at(j);
                        dof->setEqsID(i, eqs_id);
                        eqs_id++;
                    }
                }
            }
            _total_pt = eqs_id;
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

    void getListOfEqsID(const std::vector<size_t> &ele_node_ids, const std::vector<size_t> &ele_node_size_order, std::vector<size_t> &local_dofmap) const
    {
        const size_t n_dof = getNumberOfDof();
        if (n_dof<2) {
            const DofMap *dofMap = getDofMap(0);
            dofMap->getListOfEqsID(ele_node_ids, ele_node_size_order[dofMap->getOrder()-1], local_dofmap);
        } else {
            for (size_t j=0; j<n_dof; j++) {
                std::vector<size_t> tmp_map;
                const DofMap *dofMap = getDofMap(j);
                dofMap->getListOfEqsID(ele_node_ids, ele_node_size_order[dofMap->getOrder()-1], tmp_map);
                local_dofmap.insert(local_dofmap.end(), tmp_map.begin(), tmp_map.end());
            }
        }
    }
};

}
