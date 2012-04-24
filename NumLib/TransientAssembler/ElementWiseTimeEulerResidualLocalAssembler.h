
#pragma once

#include <cassert>

#include "MeshLib/Core/IElement.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "IElementWiseTransientResidualLocalAssembler.h"


namespace NumLib
{

/**
 * \brief Euler scheme element assembler for time ODE formulations
 *
 * @tparam  T_USER_ASSEMBLY 	User-given assembler
 */
template <class T_USER_ASSEMBLY>
class ElementWiseTimeEulerResidualLocalAssembler : public IElementWiseTransientResidualLocalAssembler
{
public:
	ElementWiseTimeEulerResidualLocalAssembler(T_USER_ASSEMBLY &a) : _time_ode(&a), _theta(1.0)
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
    /// @param time			time step
    /// @param e			element
    /// @param local_u_n1	guess of current time step value
    /// @param local_u_n	previous time step value
    /// @param eqs			local algebraic equation
    virtual void assembly(const TimeStep &time, MeshLib::IElement &e, const LocalVectorType &local_u_n1, const LocalVectorType &local_u_n, LocalVectorType &local_r)
    {
        const double delta_t = time.getTimeStepSize();
        const size_t n_dof = local_r.size();

        LocalMatrixType M(n_dof, n_dof);
        LocalMatrixType K(n_dof, n_dof);
        LocalVectorType F(.0, n_dof);
        M = .0;
        K = .0;

        // get M,K,F in M du/dt + K = F
        _time_ode->assembly(time, e, local_u_n1, local_u_n, M, K, F);

        //std::cout << "M="; M.write(std::cout); std::cout << std::endl;
        //std::cout << "K="; K.write(std::cout); std::cout << std::endl;

        LocalMatrixType TMP_M(n_dof, n_dof);
        LocalMatrixType TMP_M2(n_dof, n_dof);
        LocalVectorType TMP_V(.0, n_dof);

        // evaluate r: r = (1/dt M + theta K) * u1 -(1/dt M - (1-theta) K) * u0 - F
        // r = (1/dt M + theta K) * u1
        TMP_M = M;
        TMP_M *= 1.0/delta_t;
        TMP_M2 = TMP_M;
        TMP_M = K;
        TMP_M *= _theta;
        TMP_M2 += TMP_M;
        TMP_V = .0;
        TMP_M2.axpy(1.0, &((LocalVectorType)local_u_n1)[0], .0, &TMP_V[0]);
        local_r += TMP_V;
        // r -= (1/dt M - (1-theta) K) u0
        TMP_M = M;
        TMP_M *= 1.0/delta_t;
        TMP_M2 = TMP_M;
        TMP_M = K;
        TMP_M *= - (1.-_theta);
        TMP_M2 += TMP_M;
        TMP_V = .0;
        TMP_M2.axpy(1.0, &((LocalVectorType)local_u_n)[0], .0, &TMP_V[0]);
        local_r -= TMP_V;
        // r -= F
        local_r -= F;
    }

private:
    T_USER_ASSEMBLY* _time_ode;
    double _theta;
};



} //end
