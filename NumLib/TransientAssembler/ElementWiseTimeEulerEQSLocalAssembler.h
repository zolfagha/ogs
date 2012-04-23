
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
 * @tparam  T_USER_ASSEMBLY 	User-given assembler
 */
template <class T_USER_ASSEMBLY>
class ElementWiseTimeEulerEQSLocalAssembler : public IElementWiseTransientLinearEQSLocalAssembler
{
public:
    ElementWiseTimeEulerEQSLocalAssembler(T_USER_ASSEMBLY &a) : _time_ode(&a), _theta(1.0)
    {
    };

    virtual ~ElementWiseTimeEulerEQSLocalAssembler() {};

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
    virtual void assembly(const TimeStep &time, MeshLib::IElement &e, const LocalVectorType &local_u_n1, const LocalVectorType &local_u_n, LocalEquationType &eqs)
    {
        const double delta_t = time.getTimeStepSize();
        const size_t n_dof = eqs.getDimension();

        LocalMatrixType M(n_dof, n_dof);
        LocalMatrixType K(n_dof, n_dof);
        LocalVectorType F(n_dof, .0);
        M = .0;
        K = .0;

        _time_ode->assembly(time, e, local_u_n1, local_u_n, M, K, F);

        //std::cout << "M="; M.write(std::cout); std::cout << std::endl;
        //std::cout << "K="; K.write(std::cout); std::cout << std::endl;

        LocalMatrixType *localA = eqs.getA();
        double *localRHS = eqs.getRHS();
        LocalMatrixType TMP_M(n_dof, n_dof);
        LocalMatrixType TMP_M2(n_dof, n_dof);

        // A = 1/dt M + theta K
        TMP_M = M;
        TMP_M *= 1.0/delta_t;
        (*localA) = TMP_M;
        TMP_M = K;
        TMP_M *= _theta;
        (*localA) += TMP_M;
        // RHS = (1/dt M - (1-theta) K) u0 + F
        TMP_M = M;
        TMP_M *= 1.0/delta_t;
        TMP_M2 = TMP_M;
        TMP_M = K;
        TMP_M *= - (1.-_theta);
        TMP_M2 += TMP_M;
        TMP_M2.axpy(1.0, &local_u_n[0], .0, localRHS);
        for (size_t i=0; i<n_dof; i++)
            localRHS[i] += F[i];
    }

private:
    T_USER_ASSEMBLY* _time_ode;
    double _theta;
};



} //end
