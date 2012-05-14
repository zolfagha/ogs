
#pragma once

#include "DiscreteLib/Core/IDiscreteVector.h"

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
	typedef DiscreteLib::IDiscreteVector<double> GlobalVector;

	virtual ~IDiscreteLinearEquationAssembler() {};

    /// assembly
    virtual void assembly( MeshLib::IMesh &msh, DofEquationIdTable &dofManager, MathLib::ILinearEquations &eqs) = 0;
};

}
