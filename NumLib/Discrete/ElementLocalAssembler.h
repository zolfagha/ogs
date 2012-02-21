
#pragma once

#include "MathLib/LinAlg/LinearEquations/DenseLinearEquations.h"

#include "MeshLib/Core/IMesh.h"
#include "FemLib/Function/FemFunction.h"
#include "NumLib/TimeStepping/TimeStep.h"


namespace NumLib
{

/**
 * \brief Interface of all element local assembler classes
 */
class IElemenetLocalAssembler
{
public:
    /// assemble a local linear equation for the given element
    virtual void assembly(MeshLib::IElement &e, MathLib::DenseLinearEquations &eqs) = 0;
};

class ITransientElemenetLocalAssembler
{
public:
    /// assemble a local linear equation for the given element
    virtual void assembly(const TimeStep &time,  MeshLib::IElement &e, MathLib::DenseLinearEquations &eqs) = 0;
};

/**
* \brief Abstract class of element assembler for time ODE formulations with FEM
 */
class AbstractTimeODEFemElementAssembler
{
public:
    AbstractTimeODEFemElementAssembler() : _fem_func(0) {};

    void setFemFunction(FemLib::FemNodalFunctionScalar& fe)
    {
        _fem_func = &fe;
    }

    void assembly(const TimeStep &time, MeshLib::IElement &e, MathLib::DenseLinearEquations::MatrixType &M, MathLib::DenseLinearEquations::MatrixType &K,  MathLib::DenseLinearEquations::VectorType &F)
    {
        FemLib::IFiniteElement* fe = _fem_func->getFiniteElement(e);
        assemblyFE(time, *fe, M, K, F);
    }

private:
    FemLib::FemNodalFunctionScalar* _fem_func;

protected:
    virtual void assemblyFE(const TimeStep &time, FemLib::IFiniteElement &fe, MathLib::DenseLinearEquations::MatrixType &M, MathLib::DenseLinearEquations::MatrixType &K,  MathLib::DenseLinearEquations::VectorType &F) = 0;
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
        _theta = v;
    }

    virtual void assembly(const TimeStep &time, MeshLib::IElement &e, MathLib::DenseLinearEquations &eqs)
    {
        double delta_t = time.getTimeStepSize();
        std::vector<double> local_u_n;

        const size_t n_dof = eqs.getDimension();

        MathLib::DenseLinearEquations::MatrixType M(n_dof, n_dof);
        MathLib::DenseLinearEquations::MatrixType K(n_dof, n_dof);
        MathLib::DenseLinearEquations::VectorType F(n_dof);
        _time_ode.assembly(time, e, M, K, F);

        MathLib::DenseLinearEquations::MatrixType *localA = eqs.getA();
        double *localRHS = eqs.getRHS();
        MathLib::DenseLinearEquations::MatrixType TMP_M(n_dof, n_dof);
        MathLib::DenseLinearEquations::MatrixType TMP_M2(n_dof, n_dof);
        TMP_M = M;
        TMP_M *= 1.0/delta_t;
        (*localA) = TMP_M;
        TMP_M = K;
        TMP_M *= _theta;
        (*localA) += TMP_M;
        TMP_M = M;
        TMP_M *= 1.0/delta_t;
        TMP_M2 = TMP_M;
        TMP_M = K;
        TMP_M *= - (1.-_theta);
        TMP_M2 += TMP_M;
        TMP_M2.axpy(1.0, &local_u_n[0], .0, localRHS);
    }
};


} //end
