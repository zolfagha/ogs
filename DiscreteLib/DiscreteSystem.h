
#pragma once

#include "Base/CodingTools.h"
#include "Base/BidirectionalMap.h"

#include "MathLib/LinAlg/LinearEquations/ILinearEquations.h"

#include "MeshLib/Core/IMesh.h"

#include "DiscreteVector.h"
#include "DiscreteLinearEquation.h"


namespace DiscreteLib
{

/**
 * \brief Discrete system based on a single mesh
 */
class DiscreteSystem
{
public:
    /// @param msh a mesh to represent discrete systems by nodes or elements or edges
    DiscreteSystem(MeshLib::IMesh& msh) : _msh(&msh) {};

    virtual ~DiscreteSystem()
    {
        Base::releaseObjectsInStdVector(_vec_linear_sys);
        Base::releaseObjectsInStdVector(_vec_vectors);
    }

    MeshLib::IMesh* getMesh() const { return _msh; };

    /// create a new linear equation
    template<class T_LINEAR_SOLVER, class T_SPARSITY_BUILDER>
    IDiscreteLinearEquation* createLinearEquation(T_LINEAR_SOLVER &linear_solver, DofMapManager &dofManager)
    {
        _vec_linear_sys.push_back(new TemplateMeshBasedDiscreteLinearEquation<T_LINEAR_SOLVER, T_SPARSITY_BUILDER>(*_msh, linear_solver, dofManager));
        //return _vec_linear_sys.size()-1;
        return _vec_linear_sys.back();
    }

    //IDiscreteLinearEquation* getLinearEquation(size_t i)
    //{
    //    return _vec_linear_sys[i];
    //}

    template<typename T>
    DiscreteVector<T>* createVector(const size_t &n) 
    {
        DiscreteVector<T>* v = new DiscreteVector<T>(n);
        _vec_vectors.push_back(v);
        return v;
    };

private:
    DISALLOW_COPY_AND_ASSIGN(DiscreteSystem);

    // discretization
    MeshLib::IMesh* _msh;

    // linear equations
    std::vector<AbstractMeshBasedDiscreteLinearEquation*> _vec_linear_sys;

protected:
    // vector
    std::vector<IDiscreteVector*> _vec_vectors;
};

/**
 * \brief Nodally decomposed discrete system
 */
class NodeDecomposedDiscreteSystem : public DiscreteSystem
{
public:
    NodeDecomposedDiscreteSystem(MeshLib::IMesh &local_msh, Base::BidirectionalMap<size_t, size_t> &msh_node_id_mapping, std::set<size_t> &ghost_nodes)
        : DiscreteSystem(local_msh), _map_global2local_node_id(&msh_node_id_mapping), _ghost_nodes(ghost_nodes), _n_global_nodes(0)
    {
    }

    virtual ~NodeDecomposedDiscreteSystem()
    {
    }

    size_t getGlobalNumberOfNodes() const {return _n_global_nodes; };

    template<typename T>
    OMPDecomposedMasterVector<T>* createVector(const size_t &n) 
    {
        OMPDecomposedMasterVector<T>* v = new OMPDecomposedMasterVector<T>(0, n);
        _vec_vectors.push_back(v);
        return v;
    };

    template<typename T>
    OMPDecomposedMasterVector<T>* getVector(const size_t &i) 
    {
        return (OMPDecomposedMasterVector<T>*)_vec_vectors[n];
    };


private:
    DISALLOW_COPY_AND_ASSIGN(NodeDecomposedDiscreteSystem);

    // discretization
    Base::BidirectionalMap<size_t, size_t>* _map_global2local_node_id;
    std::set<size_t> _ghost_nodes;

    // 
    size_t _n_global_nodes;

};

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
    OMPDecomposedMasterVector<T>* createVector(const size_t &n) 
    {
        OMPDecomposedMasterVector<T>* v = new OMPDecomposedMasterVector<T>(0, n);
        _vec_vectors.push_back(v);
        return v;
    };

    template<typename T>
    OMPDecomposedMasterVector<T>* getVector(const size_t &i) 
    {
        return (OMPDecomposedMasterVector<T>*)_vec_vectors[n];
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
    OMPMasterNodeDecomposedDiscreteSystem(MeshLib::IMesh &global_msh, size_t n_dom)
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
    IDiscreteLinearEquation* createLinearEquation(T_LINEAR_SOLVER &linear_solver, DofMapManager &dofManager)
    {
        _vec_linear_sys.push_back(new TemplateMeshBasedDiscreteLinearEquation<T_LINEAR_SOLVER, T_SPARSITY_BUILDER>(*_msh, linear_solver, dofManager));
        //return _vec_linear_sys.size()-1;
        return _vec_linear_sys.back();
    }

    template<typename T>
    OMPDecomposedMasterVector<T>* createVector(const size_t &n) 
    {
        OMPDecomposedMasterVector<T>* v = new OMPDecomposedMasterVector<T>(0, n);
        _vec_vectors.push_back(v);
        return v;
    };

    template<typename T>
    OMPDecomposedMasterVector<T>* getVector(const size_t &i) 
    {
        return (OMPDecomposedMasterVector<T>*)_vec_vectors[n];
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
