/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file IElementWiseTransientJacobianLocalAssembler.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

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
 * \brief Interface class for element local assembler for transient case
 */
class IElementWiseTransientJacobianLocalAssembler
{
public:
    virtual ~IElementWiseTransientJacobianLocalAssembler() {};

    /// assemble a local Jacobian matrix for the given element
    /// @param time            time step
    /// @param e            element
    /// @param local_u_n1    guess of current time step value
    /// @param local_u_n    previous time step value
    /// @param local_J        local Jacobian
    virtual void assembly(const TimeStep &time,  const MeshLib::IElement &e, const LocalVector &local_u_n1, const LocalVector &local_u_n, LocalMatrix &local_J) = 0;
};


class DummyElementWiseTransientJacobianLocalAssembler : public IElementWiseTransientJacobianLocalAssembler
{
public:
    virtual ~DummyElementWiseTransientJacobianLocalAssembler() {};

    /// assemble a local Jacobian matrix for the given element
    /// @param time            time step
    /// @param e            element
    /// @param local_u_n1    guess of current time step value
    /// @param local_u_n    previous time step value
    /// @param local_J        local Jacobian
    virtual void assembly(const TimeStep &,  const MeshLib::IElement &, const LocalVector &, const LocalVector &, LocalMatrix &) {};
};


} //end
