
#pragma once

#include "NumLib/DataType.h"


namespace MeshLib
{
class IElement;
};

namespace NumLib
{

class TimeStep;

/**
 * \brief Interface class of element assembler for time ODE formulations
 */
class IElementWiseTimeODELocalAssembler
{
public:
    virtual ~IElementWiseTimeODELocalAssembler() {};

    /// assemble components of a local linear equation for the given element
    /// @param time		time step
    /// @param e		element
    /// @param local_u_n1	guess of current time step value
    /// @param local_u_n	previous time step value
    /// @param M		mass matrix
    /// @param K		laplace matrix
    /// @param F		source/sink terms
    virtual void assembly(const TimeStep &time, MeshLib::IElement &e, const LocalVector &local_u_n1, const LocalVector &local_u_n, LocalMatrix &M, LocalMatrix &K, LocalVector &F)  = 0;
};

} //end
