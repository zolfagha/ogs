
#pragma once

#include "DiscreteLib/Core/LocalDataType.h"

namespace MathLib
{
class DenseLinearEquations;
}

namespace MeshLib
{
class IElement;
}

namespace DiscreteLib
{

/**
 * \brief Interface of all element local assembler classes
 */
class IElemenetWiseLinearEquationLocalAssembler
{
public:
	virtual ~IElemenetWiseLinearEquationLocalAssembler() {};

    /// assemble a local linear equation for the given element
    virtual void assembly(MeshLib::IElement &e, LocalEquation &eqs) = 0;
};

} //end
