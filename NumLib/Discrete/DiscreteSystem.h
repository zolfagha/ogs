
#pragma once

#include "Base/CodingTools.h"

#include "MathLib/LinAlg/LinearEquations/ILinearEquations.h"

#include "MeshLib/Core/IMesh.h"

#include "NumLib/Discrete/DiscreteLinearEquation.h"


namespace NumLib
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
    }

    /// create a new linear equation
    template<class T_LINEAR_SOLVER>
    size_t addLinearEquation(T_LINEAR_SOLVER &linear_solver)
    {
        _vec_linear_sys.push_back(new TemplateMeshBasedDiscreteLinearEquation<T_LINEAR_SOLVER>(*_msh, linear_solver));
        return _vec_linear_sys.size()-1;
    }

    IDiscreteLinearEquation* getLinearEquation(size_t i)
    {
        return _vec_linear_sys[i];
    }

    //void doEachElement(int func1, int func_reduce);

    //void createField();

    MeshLib::IMesh* getMesh() const { return _msh; };

private:
    DISALLOW_COPY_AND_ASSIGN(DiscreteSystem);

    MeshLib::IMesh* _msh;
    std::vector<AbstractMeshBasedDiscreteLinearEquation*> _vec_linear_sys;

};


}