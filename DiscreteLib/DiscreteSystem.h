
#pragma once

#include "Base/CodingTools.h"
#include "Base/BidirectionalMap.h"

#include "MathLib/LinAlg/LinearEquations/ILinearEquations.h"

#include "MeshLib/Core/IMesh.h"

#include "DiscreteVector.h"
#include "DiscreteLinearEquation.h"

namespace DiscreteLib
{

class IDiscreteSystem
{
};

/**
 * \brief Discrete system based on a single mesh
 */
class DiscreteSystem : public IDiscreteSystem
{
public:
    /// @param msh a mesh to represent discrete systems by nodes or elements or edges
    DiscreteSystem(MeshLib::IMesh& msh) : _msh(&msh) {};

    virtual ~DiscreteSystem()
    {
        Base::releaseObjectsInStdVector(_vec_linear_sys);
        Base::releaseObjectsInStdVector(_vec_vectors);
    }

    /// get this mesh
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

} //end
