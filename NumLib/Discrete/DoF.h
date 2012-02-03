
#pragma once

#include <vector>
#include <map>
#include "Base/MemoryTools.h"
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

public:
    DofMap(size_t discrete_points_size) : _map_node_id2eqs_id(discrete_points_size)
    {
    }

    void getListOfEqsID( const std::vector<size_t>& vec_pt_id, std::vector<size_t>& vec_eqs_id ) const
    {
        vec_eqs_id.resize(vec_pt_id.size());
        for (size_t i=0; i<vec_pt_id.size(); i++) {
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
        Base::destroyStdVectorWithPointers(_map_var2dof);
    }

    size_t addDoF(size_t discrete_points_size, size_t mesh_id=0)
    {
        DofMap *dof = new DofMap(discrete_points_size);
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
};

}
