
#pragma once

#include <vector>

#include "MathLib/LinAlg/LinearEquations/DenseLinearEquations.h"

#include "MeshLib/Core/IMesh.h"


namespace DiscreteLib
{

/**
 * \brief Interface of all element local assembler classes
 */
class IElemenetLocalAssembler
{
public:
	virtual ~IElemenetLocalAssembler() {};
    /// assemble a local linear equation for the given element
    virtual void assembly(MeshLib::IElement &e, MathLib::DenseLinearEquations &eqs) = 0;
};

} //end
