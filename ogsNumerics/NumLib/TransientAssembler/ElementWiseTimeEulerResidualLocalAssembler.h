/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ElementWiseTimeEulerResidualLocalAssembler.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <cassert>

#include "logog.hpp"


#include "MeshLib/Core/IElement.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "IElementWiseTransientResidualLocalAssembler.h"


namespace NumLib
{

/**
 * \brief Euler scheme element assembler for time ODE formulations
 *
 * @tparam  T_USER_ASSEMBLY     User-given assembler
 */
class ElementWiseTimeEulerResidualLocalAssembler : public IElementWiseTransientResidualLocalAssembler
{
public:
    ElementWiseTimeEulerResidualLocalAssembler() : _theta(1.0)
    {
    };

    virtual ~ElementWiseTimeEulerResidualLocalAssembler() {};

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
    virtual void assembly(const TimeStep &time, const MeshLib::IElement &e, const LocalVector &local_u_n1, const LocalVector &local_u_n, LocalVector &local_r)
    {
        const double delta_t = time.getTimeStepSize();
        const size_t n_dof = local_r.size();

        LocalMatrix M(n_dof, n_dof);
        LocalMatrix K(n_dof, n_dof);
        LocalVector F(n_dof);
        M *= .0;
        K *= .0;
        F *= .0;

        // get M,K,F in M du/dt + K = F
        assembleODE(time, e, local_u_n1, local_u_n, M, K, F);

        LocalMatrix TMP_M(n_dof, n_dof);
        LocalVector TMP_V(n_dof);

        // evaluate r: r = (1/dt M + theta K) * u1 -(1/dt M - (1-theta) K) * u0 - F
        // r = (1/dt M + theta K) * u1
        TMP_M = 1.0/delta_t * M + _theta * K;
        local_r = TMP_M * local_u_n1;
        // r -= (1/dt M - (1-theta) K) u0
        TMP_M = 1.0/delta_t * M - (1.-_theta) * K;
        local_r.noalias() -= TMP_M * local_u_n;
        // r -= F
        local_r.noalias() -= F;

#if 0
        if (e.getID()<1) {
            std::cout << "Element: " << e.getID() << std::endl;
            std::cout << "M=" << M << std::endl;
            std::cout << "K=" << K << std::endl;
            std::cout << "u1=" << local_u_n1 << std::endl;
            std::cout << "u0=" << local_u_n << std::endl;
            TMP_M = 1.0/delta_t * M + _theta * K;
            std::cout << "A=" << TMP_M << std::endl;
            TMP_M = 1.0/delta_t * M - (1.-_theta) * K;
            std::cout << "B=" << TMP_M << std::endl;
            std::cout << "r=" << local_r << std::endl;
        }

#endif
    }

protected:
    virtual void assembleODE(const TimeStep &time, const MeshLib::IElement &e, const LocalVector &local_u_n1, const LocalVector &local_u_n, LocalMatrix &M, LocalMatrix &K, LocalVector &F)  = 0;

private:
    double _theta;
};



} //end
