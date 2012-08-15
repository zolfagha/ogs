/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file OMPLocalNodeDecomposedDiscreteSystem.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <set>

#include "BaseLib/CodingTools.h"
#include "BaseLib/BidirectionalMap.h"
#include "DiscreteLib/Core/DiscreteSystem.h"
#include "OMPGlobalDiscreteVector.h"

namespace DiscreteLib
{

/**
 * \brief Local discrete system for OpenMP node-based decomposition
 */
class OmpNodeDdcLocalDiscreteSystem : public DiscreteSystem
{
public:
    /// constructor
    /// \param local_msh            a local mesh object
    /// \param msh_node_id_mapping  mapping table from global node id to local node id
    /// \param ghost_nodes          a list of ghost nodes
    OmpNodeDdcLocalDiscreteSystem(MeshLib::IMesh &local_msh, BaseLib::BidirectionalMap<size_t, size_t> &msh_node_id_mapping, std::set<size_t> &ghost_nodes)
        : DiscreteSystem(local_msh), _map_global2local_node_id(&msh_node_id_mapping), _ghost_nodes(ghost_nodes), _n_global_nodes(0)
    {
    }

    ///
    virtual ~OmpNodeDdcLocalDiscreteSystem() { }

    /// return the number of all nodes
    size_t getGlobalNumberOfNodes() const {return _n_global_nodes; };

    /// create a discrete vector
    template<typename T>
    OMPGlobalDiscreteVector<T>* createVector(const size_t &n) 
    {
        OMPGlobalDiscreteVector<T>* v = new OMPGlobalDiscreteVector<T>(0, n);
        DiscreteSystem::_data.addVector(v);
        return v;
    };

    /// get a stored discrete vector
    template<typename T>
    OMPGlobalDiscreteVector<T>* getVector(const size_t &i) 
    {
        return (OMPGlobalDiscreteVector<T>*)DiscreteSystem::_data.getVector(i);
    };


private:
    DISALLOW_COPY_AND_ASSIGN(OmpNodeDdcLocalDiscreteSystem);

private:
    BaseLib::BidirectionalMap<size_t, size_t>* _map_global2local_node_id;
    std::set<size_t> _ghost_nodes;
    size_t _n_global_nodes;

};

} //end
