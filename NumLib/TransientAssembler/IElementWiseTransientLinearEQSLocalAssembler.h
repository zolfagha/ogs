
#pragma once

#include "MathLib/LinAlg/LinearEquations/DenseLinearEquations.h"
#include "NumLib/DataType.h"


namespace MeshLib
{
class IElement;
};

namespace NumLib
{

class TimeStep;

/**
 * \brief Interface class for element local assembler for transient case
 */
class IElementWiseTransientLinearEQSLocalAssembler
{
public:
	typedef MathLib::DenseLinearEquations LocalEquationType;

    virtual ~IElementWiseTransientLinearEQSLocalAssembler() {};

    /// assemble a local linear equation for the given element
    /// @param time			time step
    /// @param e			element
    /// @param local_u_n1	guess of current time step value
    /// @param local_u_n	previous time step value
    /// @param eqs			local algebraic equation
    virtual void assembly(const TimeStep &time,  MeshLib::IElement &e, const LocalVector &local_u_n1, const LocalVector &local_u_n, LocalEquationType &eqs) = 0;
};

} //end
