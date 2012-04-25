
#pragma once

#include "DiscreteLib/Core/DiscreteSystem.h"
#include "DiscreteLib/LinearEquation/MeshBasedDiscreteLinearEquation.h"
#include "OMPDiscreteVector.h"

namespace DiscreteLib
{

class OMPLocalNodeDecomposedDiscreteSystem : public DiscreteSystem
{
public:
    OMPLocalNodeDecomposedDiscreteSystem(MeshLib::IMesh &local_msh, Base::BidirectionalMap<size_t, size_t> &msh_node_id_mapping, std::set<size_t> &ghost_nodes)
        : DiscreteSystem(local_msh), _map_global2local_node_id(&msh_node_id_mapping), _ghost_nodes(ghost_nodes), _n_global_nodes(0)
    {
    }

    virtual ~OMPLocalNodeDecomposedDiscreteSystem()
    {
    }

    size_t getGlobalNumberOfNodes() const {return _n_global_nodes; };

    template<typename T>
    OMPGlobalDiscreteVector<T>* createVector(const size_t &n) 
    {
        OMPGlobalDiscreteVector<T>* v = new OMPGlobalDiscreteVector<T>(0, n);
	DiscreteSystem::_data.addVector(v);
        return v;
    };

    template<typename T>
    OMPGlobalDiscreteVector<T>* getVector(const size_t &i) 
    {
    	return (OMPGlobalDiscreteVector<T>*)DiscreteSystem::_data.getLinearEquation(i);
    };


private:
    DISALLOW_COPY_AND_ASSIGN(OMPLocalNodeDecomposedDiscreteSystem);

    // discretization
    Base::BidirectionalMap<size_t, size_t>* _map_global2local_node_id;
    std::set<size_t> _ghost_nodes;

    // 
    size_t _n_global_nodes;

};

/**
 * \brief Nodally decomposed discrete system
 */
class OMPMasterNodeDecomposedDiscreteSystem : public DiscreteSystem
{
public:
    OMPMasterNodeDecomposedDiscreteSystem(MeshLib::IMesh &global_msh, size_t /*n_dom*/)
        : DiscreteSystem(global_msh)
    {
    }

    virtual ~OMPMasterNodeDecomposedDiscreteSystem()
    {
    }

    size_t getGlobalNumberOfNodes() const {return _n_global_nodes; };

    OMPLocalNodeDecomposedDiscreteSystem* createLocal(MeshLib::IMesh &local_msh, Base::BidirectionalMap<size_t, size_t> &msh_node_id_mapping, std::set<size_t> &ghost_nodes)
    {
        OMPLocalNodeDecomposedDiscreteSystem* local = new OMPLocalNodeDecomposedDiscreteSystem(local_msh, msh_node_id_mapping, ghost_nodes);
        _vec_local_dis.push_back(local);
        return local;
    }

    OMPLocalNodeDecomposedDiscreteSystem* getLocal(size_t i) {return _vec_local_dis[i];};
    size_t getNumberOfLocal() const {return _vec_local_dis.size();};

    /// create a new linear equation
    template<class T_LINEAR_SOLVER, class T_SPARSITY_BUILDER>
    IDiscreteLinearEquation* createLinearEquation(T_LINEAR_SOLVER &linear_solver, DofEquationIdTable &dofManager)
    {
	
	TemplateMeshBasedDiscreteLinearEquation<T_LINEAR_SOLVER, T_SPARSITY_BUILDER>* eq = new TemplateMeshBasedDiscreteLinearEquation<T_LINEAR_SOLVER, T_SPARSITY_BUILDER>(*_msh, linear_solver, dofManager);
	DiscreteSystem::_data.addLinearEquation(eq);
        //return _vec_linear_sys.size()-1;
        return eq;
    }

    template<typename T>
    OMPGlobalDiscreteVector<T>* createVector(const size_t &n) 
    {
        OMPGlobalDiscreteVector<T>* v = new OMPGlobalDiscreteVector<T>(0, n);
        DiscreteSystem::_data.addVector(v);
        return v;
    };

    template<typename T>
    OMPGlobalDiscreteVector<T>* getVector(const size_t &i) 
    {
        return (OMPGlobalDiscreteVector<T>*)DiscreteSystem::_data.getVector(i);
    };


private:
    DISALLOW_COPY_AND_ASSIGN(OMPMasterNodeDecomposedDiscreteSystem);

    std::vector<OMPLocalNodeDecomposedDiscreteSystem*> _vec_local_dis;
    // discretization
    Base::BidirectionalMap<size_t, size_t>* _map_global2local_node_id;
    std::set<size_t> _ghost_nodes;

    // 
    size_t _n_global_nodes;

};


} //end
