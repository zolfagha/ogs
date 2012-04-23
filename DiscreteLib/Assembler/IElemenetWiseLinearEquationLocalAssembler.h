
#pragma once

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
    virtual void assembly(MeshLib::IElement &e, MathLib::DenseLinearEquations &eqs) = 0;
};

} //end
