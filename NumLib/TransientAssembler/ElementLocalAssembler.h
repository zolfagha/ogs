
#pragma once

#include <vector>

#include "MathLib/LinAlg/LinearEquations/DenseLinearEquations.h"

#include "NumLib/TimeStepping/TimeStep.h"


namespace MeshLib
{
class IElement;
};

namespace NumLib
{

/**
 * \brief Interface class for element local assembler for transient case
 */
class ITransientElemenetLocalAssembler
{
public:
	typedef MathLib::DenseLinearEquations LocalEquationType;
	typedef LocalEquationType::VectorType LocalVectorType;

    virtual ~ITransientElemenetLocalAssembler() {};

    /// assemble a local linear equation for the given element
    /// @param time			time step
    /// @param e			element
    /// @param local_u_n1	guess of current time step value
    /// @param local_u_n	previous time step value
    /// @param eqs			local algebraic equation
    virtual void assembly(const TimeStep &time,  MeshLib::IElement &e, const LocalVectorType &local_u_n1, const LocalVectorType &local_u_n, LocalEquationType &eqs) = 0;
};

/**
 * \brief Interface class of element assembler for time ODE formulations
 */
class ITimeODEElementAssembler
{
public:
	typedef MathLib::DenseLinearEquations LocalEquationType;
	typedef LocalEquationType::MatrixType LocalMatrixType;
	typedef LocalEquationType::VectorType LocalVectorType;

    virtual ~ITimeODEElementAssembler() {};

    /// assemble components of a local linear equation for the given element
    /// @param time		time step
    /// @param e		element
    /// @param local_u_n1	guess of current time step value
    /// @param local_u_n	previous time step value
    /// @param M		mass matrix
    /// @param K		laplace matrix
    /// @param F		source/sink terms
    virtual void assembly(const TimeStep &time, MeshLib::IElement &e, const LocalVectorType &local_u_n1, const LocalVectorType &local_u_n, LocalMatrixType &M, LocalMatrixType &K, LocalVectorType &F)  = 0;
};

} //end
