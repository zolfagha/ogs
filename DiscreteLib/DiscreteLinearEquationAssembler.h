
#pragma once

#include "ElementLocalAssembler.h"
#include "DoF.h"

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
 * \brief Interface of discrete system assembler classes
 */
class IDiscreteLinearEquationAssembler
{
public:
    /// assembly
    virtual void assembly( MeshLib::IMesh &msh, DofMapManager &dofManager, MathLib::ILinearEquations &eqs) = 0;
};

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
    void assembly(MeshLib::IMesh &msh, DofMapManager &dofManager, MathLib::ILinearEquations &eqs);

private:
    IElemenetLocalAssembler* _e_assembler;
};


class NodeBasedAssembler  : public IDiscreteLinearEquationAssembler
{
public:
    void assembly( MeshLib::IMesh& _msh, MathLib::ILinearEquations& _eqs) {};
};

class EdgeBasedAssembler  : public IDiscreteLinearEquationAssembler
{
public:
    void assembly( MeshLib::IMesh& _msh, MathLib::ILinearEquations& _eqs) {};
};

}
