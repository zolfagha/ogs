
#pragma once

#include "MathLib/Function/Function.h"
#include "MathLib/LinAlg/MatrixBase.h"
#include "Convergence.h"

namespace MathLib
{
class NewtoRaphson
{
public:
	NewtoRaphson()
    {
        _max_itr_count = 100;
    }

    virtual ~NewtoRaphson()
    {
    }

    template <class F_RESIDUAL, class F_JACOBIAN, class T_MATRIX, class T_VECTOR, class T_LINEAR_SOLVER, class T_CONVERGENCE>
    int solve(F_RESIDUAL &funcResidual, F_JACOBIAN &funcJacobian, T_VECTOR &x0, T_VECTOR &x_new)
    {
    	T_CONVERGENCE convergence_check;
    	T_MATRIX jacobian;
    	T_VECTOR r;
    	T_VECTOR x_prev = x0;
    	T_VECTOR dx;
        T_LINEAR_SOLVER linear_solver;

    	for (size_t i=0; i<_max_itr_count; i++) {
        	funcJacobian.eval(x_prev, jacobian);
            funcResidual.eval(x_prev, r);
            linear_solver.solve(jacobian, r, dx);
            x_new = x_prev + dx;
            if (convergence_check.isConverged(dx))
                return 0;
        }

        return -1; //not converged
    }

private:
    size_t _max_itr_count;
};


template<class T_LINEAR_EQUATION>
class TemplateNewtoRaphson
{
public:
	TemplateNewtoRaphson(T_LINEAR_EQUATION &linear_eqs)
    {
        _max_itr_count = 100;
        _linear_eqs = &linear_eqs;
    }

    virtual ~TemplateNewtoRaphson()
    {
    }

    template <class F_RESIDUAL, class F_JACOBIAN, class T_MATRIX, class T_VECTOR, class T_LINEAR_SOLVER, class T_CONVERGENCE>
    int solve(F_RESIDUAL &funcResidual, F_JACOBIAN &funcJacobian, T_VECTOR &x0, T_VECTOR &x_new)
    {
    	T_CONVERGENCE convergence_check;
    	T_MATRIX jacobian;
    	T_VECTOR r;
    	T_VECTOR x_prev = x0;
    	T_VECTOR dx;
        T_LINEAR_SOLVER linear_solver;

    	for (size_t i=0; i<_max_itr_count; i++) {
        	funcJacobian.eval(x_prev, jacobian);
            funcResidual.eval(x_prev, r);
            linear_solver.solve(jacobian, r, dx);
            x_new = x_prev + dx;
            if (convergence_check.isConverged(dx))
                return 0;
        }

        return -1; //not converged
    }

private:
    size_t _max_itr_count;
    T_LINEAR_EQUATION* _linear_eqs;

};

}
