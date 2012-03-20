
#pragma once

#include "DiscreteLinearEquationAssembler.h"

namespace MeshLib 
{
class IMesh;
}

namespace MathLib
{
    class ILinearEquations;
}

namespace DiscreteLib
{

class IElemenetLocalAssembler;

/**
 * \brief Element-based discrete system assembler classes
 */
class ElementBasedAssembler : public IDiscreteLinearEquationAssembler
{
public:
    ///
    ElementBasedAssembler(IElemenetLocalAssembler &a) : _e_assembler(&a) {};

    /// Conduct the element by element assembly procedure
    ///
    /// @param msh Mesh
    /// @param dofManager Dof map manager
    /// @param list_dofId List of Dof IDs used in this problem
    /// @param eqs Linear equation solver
    void assembly(MeshLib::IMesh &msh, DofEquationIdTable &dofManager, MathLib::ILinearEquations &eqs);

private:
    IElemenetLocalAssembler* _e_assembler;
};


}
