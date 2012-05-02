
#pragma once

#include "MathLib/LinAlg/LinearEquations/DenseLinearEquations.h"


namespace MeshLib
{
class IStencil;
};

namespace NumLib
{

class TimeStep;

/**
 * \brief Interface class for element local assembler for transient case
 */
class IStencilWiseTransientResidualLocalAssembler
{
public:
	typedef MathLib::DenseLinearEquations LocalEquationType;
	typedef LocalEquationType::MatrixType LocalMatrixType;
	typedef LocalEquationType::VectorType LocalVectorType;

    virtual ~IStencilWiseTransientResidualLocalAssembler() {};

    /// assemble a local residual for the given element
    /// @param time			time step
    /// @param e			element
    /// @param local_u_n1	guess of current time step value
    /// @param local_u_n	previous time step value
    /// @param local_r		local residual
    virtual void assembly(const TimeStep &time,  MeshLib::IStencil &e, const LocalVectorType &local_u_n1, const LocalVectorType &local_u_n, LocalVectorType &local_r) = 0;
};

} //end
