
#pragma once

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
class DofEquationIdTable;

/**
 * \brief Interface of discrete system assembler classes
 */
class IDiscreteLinearEquationAssembler
{
public:
	virtual ~IDiscreteLinearEquationAssembler() {};

    /// assembly
    virtual void assembly( MeshLib::IMesh &msh, DofEquationIdTable &dofManager, MathLib::ILinearEquations &eqs) = 0;
};

}
