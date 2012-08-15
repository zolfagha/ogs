/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file OMPMasterNodeDecomposedDiscreteSystem.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <set>
#include <vector>

#include "BaseLib/CodingTools.h"
#include "BaseLib/BidirectionalMap.h"

#include "MeshLib/Core/IMesh.h"

#include "DiscreteLib/Core/DiscreteSystem.h"
#include "DiscreteLib/LinearEquation/MeshBasedDiscreteLinearEquation.h"
#include "OMPGlobalDiscreteVector.h"
#include "OMPLocalNodeDecomposedDiscreteSystem.h"

namespace DiscreteLib
{

/**
 * \brief Discrete system for OpenMP node-based decomposition
 */
class OmpNodeDdcDiscreteSystem : public DiscreteSystem
{
public:
    /// Constructor
    ///
    /// \param global_msh   a global mesh object
    /// \param n_dom        the number of sub domains
    OmpNodeDdcDiscreteSystem(MeshLib::IMesh &global_msh, size_t /*n_dom*/)
        : DiscreteSystem(global_msh), _n_global_nodes(0)
    {
    }

    ///
    virtual ~OmpNodeDdcDiscreteSystem()
    {
        BaseLib::releaseObjectsInStdVector(_vec_local_sys);
    }

    /// get the number of global nodes
    size_t getGlobalNumberOfNodes() const {return _n_global_nodes; };

    /// create a local discrete system
    ///
    /// \param local_msh            a local mesh object
    /// \param msh_node_id_mapping  mapping table from global node id to local node id
    /// \param ghost_nodes          a list of ghost nodes
    OmpNodeDdcLocalDiscreteSystem* createLocal(MeshLib::IMesh &local_msh, BaseLib::BidirectionalMap<size_t, size_t> &msh_node_id_mapping, std::set<size_t> &ghost_nodes)
    {
        OmpNodeDdcLocalDiscreteSystem* local = new OmpNodeDdcLocalDiscreteSystem(local_msh, msh_node_id_mapping, ghost_nodes);
        _vec_local_sys.push_back(local);
        return local;
    }

    /// return a stored local discrete system
    OmpNodeDdcLocalDiscreteSystem* getLocal(size_t i) const
    {
        assert (i<_vec_local_sys.size());
        return _vec_local_sys[i];
    };

    /// return the number of local systems
    size_t getNumberOfLocal() const {return _vec_local_sys.size();};

    /// create a new linear equation
    template<class T_LINEAR_SOLVER, class T_SPARSITY_BUILDER>
    IDiscreteLinearEquation* createLinearEquation(T_LINEAR_SOLVER &linear_solver, DofEquationIdTable &dofManager)
    {
        TemplateMeshBasedDiscreteLinearEquation<T_LINEAR_SOLVER, T_SPARSITY_BUILDER>* eq;
        eq = new TemplateMeshBasedDiscreteLinearEquation<T_LINEAR_SOLVER, T_SPARSITY_BUILDER>(*_msh, linear_solver, dofManager);
        DiscreteSystem::_data.addLinearEquation(eq);
        return eq;
    }

    /// create a discrete vector
    template<typename T>
    OMPGlobalDiscreteVector<T>* createVector(const size_t &n) 
    {
        OMPGlobalDiscreteVector<T>* v = new OMPGlobalDiscreteVector<T>(0, n);
        DiscreteSystem::_data.addVector(v);
        return v;
    };

    /// return a discrete vector
    template<typename T>
    OMPGlobalDiscreteVector<T>* getVector(const size_t &i) 
    {
        assert (i<DiscreteSystem::_data.getNumberOfVectors());
        return (OMPGlobalDiscreteVector<T>*)DiscreteSystem::_data.getVector(i);
    };

private:
    DISALLOW_COPY_AND_ASSIGN(OmpNodeDdcDiscreteSystem);


private:
    std::vector<OmpNodeDdcLocalDiscreteSystem*> _vec_local_sys;
    std::set<size_t> _ghost_nodes;
    size_t _n_global_nodes;

};


} //end
