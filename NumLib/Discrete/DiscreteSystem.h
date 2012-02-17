
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
    /// @param linear_eqs Linear equation objects used to solve linear equations. 
    DiscreteSystem(MeshLib::IMesh& msh, MathLib::ILinearEquations& linear_eqs) 
    {
        _linear_sys = new MeshBasedDiscreteLinearEquation(msh, linear_eqs);
    };

    virtual ~DiscreteSystem()
    {
        if (_linear_sys) delete _linear_sys;
    }


    IDiscreteLinearEquation* getLinearEquation()
    {
        return _linear_sys;
    }

    void doEachElement(int func1, int func_reduce);

    void createField();

private:
    DISALLOW_COPY_AND_ASSIGN(DiscreteSystem);

    MeshBasedDiscreteLinearEquation* _linear_sys;

};

class MultipleMeshDiscreteSystem;


}