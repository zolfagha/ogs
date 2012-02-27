
#pragma once

#include <vector>

#include "MathLib/LinAlg/LinearEquations/DenseLinearEquations.h"

#include "MeshLib/Core/IMesh.h"
#include "NumLib/TimeStepping/TimeStep.h"


namespace NumLib
{

class ITransientElemenetLocalAssembler
{
public:
    /// assemble a local linear equation for the given element
    virtual void assembly(const TimeStep &time,  MeshLib::IElement &e, const std::vector<double> &local_u_n, MathLib::DenseLinearEquations &eqs) = 0;
};

/**
* \brief Abstract class of element assembler for time ODE formulations with FEM
 */
class ITimeODEElementAssembler
{
public:
    virtual void assembly(const TimeStep &time, MeshLib::IElement &e, MathLib::DenseLinearEquations::MatrixType &M, MathLib::DenseLinearEquations::MatrixType &K,  MathLib::DenseLinearEquations::VectorType &F)  = 0;
};

/**
* \brief Euler scheme element assembler for time ODE formulations
 */
template <class T_USER_ASSEMBLY>
class TimeEulerElementAssembler : public ITransientElemenetLocalAssembler
{
private:
    T_USER_ASSEMBLY _time_ode;
    double _theta;
public:
    TimeEulerElementAssembler(T_USER_ASSEMBLY &a) : _theta(1.0), _time_ode(a)
    {
    };

    void setTheta(double v)
    {
        assert(v>=.0 && v<=1.0);
        _theta = v;
    }

    //T_USER_ASSEMBLY* getUserAssembler() {return &_time_ode;};

    virtual void assembly(const TimeStep &time, MeshLib::IElement &e, const std::vector<double> &local_u_n, MathLib::DenseLinearEquations &eqs)
    {
        const double delta_t = time.getTimeStepSize();
        const size_t n_dof = eqs.getDimension();

        MathLib::DenseLinearEquations::MatrixType M(n_dof, n_dof);
        MathLib::DenseLinearEquations::MatrixType K(n_dof, n_dof);
        MathLib::DenseLinearEquations::VectorType F(n_dof, .0);
        M = .0;
        K = .0;

        _time_ode.assembly(time, e, M, K, F);

        MathLib::DenseLinearEquations::MatrixType *localA = eqs.getA();
        double *localRHS = eqs.getRHS();
        MathLib::DenseLinearEquations::MatrixType TMP_M(n_dof, n_dof);
        MathLib::DenseLinearEquations::MatrixType TMP_M2(n_dof, n_dof);

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
};


} //end
