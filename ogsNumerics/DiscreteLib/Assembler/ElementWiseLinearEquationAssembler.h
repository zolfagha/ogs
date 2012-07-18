
#pragma once

#include "IDiscreteLinearEquationAssembler.h"
#include "IElemenetWiseLinearEquationLocalAssembler.h"

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

/**
 * \brief Element-based discrete system assembler classes
 */
class ElementWiseLinearEquationAssembler : public IDiscreteLinearEquationAssembler
{
public:
    ///
    explicit ElementWiseLinearEquationAssembler(IElemenetWiseLinearEquationLocalAssembler &a) : _e_assembler(&a) {};

    /// Conduct the element by element assembly procedure
    ///
    /// @param msh Mesh
    /// @param dofManager Dof map manager
    /// @param list_dofId List of Dof IDs used in this problem
    /// @param eqs Linear equation solver
    void assembly(MeshLib::IMesh &msh, DofEquationIdTable &dofManager, MathLib::ILinearEquations &eqs);

private:
    IElemenetWiseLinearEquationLocalAssembler* _e_assembler;
};


}
