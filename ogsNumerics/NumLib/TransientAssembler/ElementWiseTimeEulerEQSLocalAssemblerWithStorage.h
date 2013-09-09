/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ElementWiseTimeEulerEQSLocalAssemblerWithStorage.h
 *
 * Created on 2013-09-09 by Haibing Shao
 */

#pragma once

#include <cassert>

#include "MeshLib/Core/IElement.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "IElementWiseTransientLinearEQSLocalAssembler.h"


namespace NumLib
{

/**
 * \brief Euler scheme element assembler for time ODE formulations
 *
 * @tparam  T_USER_ASSEMBLY     User-given assembler
 */
class ElementWiseTimeEulerEQSLocalAssemblerWithStorage : public IElementWiseTransientLinearEQSLocalAssembler
{
public:
    ElementWiseTimeEulerEQSLocalAssemblerWithStorage() : _theta(1.0), _stored_t(0.0)
    {
    };

    virtual ~ElementWiseTimeEulerEQSLocalAssemblerWithStorage() {};

    ///
    void setTheta(double v)
    {
        assert(v>=.0 && v<=1.0);
        _theta = v;
    }

    /// assemble a local linear equation for the given element
    /// @param time            time step
    /// @param e            element
    /// @param local_u_n1    guess of current time step value
    /// @param local_u_n    previous time step value
    /// @param eqs            local algebraic equation
    virtual void assembly(const TimeStep &time, const MeshLib::IElement &e, const DiscreteLib::DofEquationIdTable &localDofManager, const MathLib::LocalVector &local_u_n1, const MathLib::LocalVector &local_u_n, MathLib::LocalEquation &eqs);

protected:
    virtual void assembleODE(const TimeStep &time, const MeshLib::IElement &e, const MathLib::LocalVector &local_u_n1, const MathLib::LocalVector &local_u_n, MathLib::LocalMatrix &M, MathLib::LocalMatrix &K, MathLib::LocalVector &F)  = 0;

private:
    double _theta; 

    /**
      * current time value
      */
    double _stored_t; 

    /**
      * Local mass matrix
      */
    MathLib::LocalMatrix _M; 

    /**
      * Local conductance matrix
      */
    MathLib::LocalMatrix _K; 

    /**
      * Local force vector
      */
    MathLib::LocalVector _F; 
};



} //end
