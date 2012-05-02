
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
class IStencilWiseTransientJacobianLocalAssembler
{
public:
	typedef MathLib::DenseLinearEquations LocalEquationType;
	typedef LocalEquationType::MatrixType LocalMatrixType;
	typedef LocalEquationType::VectorType LocalVectorType;

    virtual ~IStencilWiseTransientJacobianLocalAssembler() {};

    /// assemble a local Jacobian matrix for the given element
    /// @param time			time step
    /// @param e			element
    /// @param local_u_n1	guess of current time step value
    /// @param local_u_n	previous time step value
    /// @param local_J		local Jacobian
    virtual void assembly(const TimeStep &time,  MeshLib::IStencil &e, const LocalVectorType &local_u_n1, const LocalVectorType &local_u_n, LocalMatrixType &local_J) = 0;
};


class DummyStencilWiseTransientJacobianLocalAssembler : public IStencilWiseTransientJacobianLocalAssembler
{
public:
	typedef MathLib::DenseLinearEquations LocalEquationType;
	typedef LocalEquationType::MatrixType LocalMatrixType;
	typedef LocalEquationType::VectorType LocalVectorType;

    virtual ~DummyStencilWiseTransientJacobianLocalAssembler() {};

    /// assemble a local Jacobian matrix for the given element
    /// @param time			time step
    /// @param e			element
    /// @param local_u_n1	guess of current time step value
    /// @param local_u_n	previous time step value
    /// @param local_J		local Jacobian
    virtual void assembly(const TimeStep &,  MeshLib::IStencil &, const LocalVectorType &, const LocalVectorType &, LocalMatrixType &) {};
};


} //end
