
#pragma once

#include "MathLib/LinAlg/LinearEquations/DenseLinearEquations.h"

#include "MeshLib/Core/IMesh.h"
#include "FemLib/Function/FemFunction.h"
#include "NumLib/TimeStepping/TimeStep.h"


namespace NumLib
{

class IElemenetLocalAssembler
{
public:
    virtual void assembly( MeshLib::IElement &e, MathLib::DenseLinearEquations &eqs) = 0;
};

class ITimeODEFemElementAssembler
{
public:
    virtual void assembly(TimeStep &time, FemLib::IFiniteElement &fe, MathLib::DenseLinearEquations::MatrixType &M, MathLib::DenseLinearEquations::MatrixType &K,  MathLib::DenseLinearEquations::VectorType &F) = 0;
};

template <class T_USER_ASSEMBLY>
class TimeODEElementAssembler
{
private:
    FemLib::FemNodalFunctionScalar* _fem_func;

public:
    TimeODEElementAssembler() {};

    TimeODEElementAssembler(FemLib::FemNodalFunctionScalar &func) : _fem_func(&func) {};

    void assembly(TimeStep &time, MeshLib::IElement &e, MathLib::DenseLinearEquations::MatrixType &M, MathLib::DenseLinearEquations::MatrixType &K,  MathLib::DenseLinearEquations::VectorType &F)
    {
        FemLib::IFiniteElement* fe = _fem_func->getFiniteElement(&e);
        T_USER_ASSEMBLY.assembly(time, *fe, M, K, F);
    }

};

template <class T_USER_ASSEMBLY>
class TimeEulerElementAssembler : public IElemenetLocalAssembler
{
private:
    TimeODEElementAssembler<T_USER_ASSEMBLY> _time_ode;
    TimeStep *_time_step;
    double _theta;
public:
    TimeEulerElementAssembler(TimeStep &time, double theta, TimeODEElementAssembler<T_USER_ASSEMBLY> &a) : _time_step(&time), _theta(theta), _time_ode(a) {};

    virtual void assembly( MeshLib::IElement &e, MathLib::DenseLinearEquations &eqs)
    {
        double delta_t = _time_step->getTimeStep();
        double theta;
        std::vector<double> local_u_n;

        const size_t n_dof = eqs.getDimension();

        MathLib::DenseLinearEquations::MatrixType M(n_dof, n_dof);
        MathLib::DenseLinearEquations::MatrixType K(n_dof, n_dof);
        MathLib::DenseLinearEquations::VectorType F(n_dof);
        _time_ode.assembly(*_time_step, e, M, K, F);

        MathLib::DenseLinearEquations::MatrixType *localA = eqs.getA();
        double *localRHS = eqs.getRHS();
        MathLib::DenseLinearEquations::MatrixType TMP_M(n_dof, n_dof);
        MathLib::DenseLinearEquations::MatrixType TMP_M2(n_dof, n_dof);
        TMP_M = M;
        TMP_M *= 1.0/delta_t;
        (*localA) = TMP_M;
        TMP_M = K;
        TMP_M *= theta;
        (*localA) += TMP_M;
        TMP_M = M;
        TMP_M *= 1.0/delta_t;
        TMP_M2 = TMP_M;
        TMP_M = K;
        TMP_M *= - (1.-theta);
        TMP_M2 += TMP_M;
        TMP_M2.axpy(1.0, &local_u_n[0], .0, localRHS);
    }
};


} //end
