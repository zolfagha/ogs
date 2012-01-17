
#pragma once

#include "Nonlinear/Convergence.h"

namespace MathLib
{

class NewtoRaphson
{
public:
    NewtoRaphson() {
        _max_itr_count = 100;
        _convergence = new ConvergenceChecker();
    }

    virtual ~NewtoRaphson() {
        if (_convergence)
            delete _convergence;
    }

    int solve(Function *funcResidual, Function *funcJacobian, size_t x_length, double* x0, double *x_new) {

        double *x_prev = x0;
        double *x_new;
        for (size_t i=0; i<_max_itr_count; i++) {
            double *jacobian = funcJacobian->eval(x_prev);
            double *r = funcResidual->eval(x_prev);
            double *dx; //solve J dx = -r
            x_new = x_prev + dx;
            if (_convergence->isConverged(dx)) 
                return 0;
        }

        return -1; //not converged
    }

    int solve(Function *function, size_t x_length, double* x0, double *x_new, IConvergenceChecker* conv_fun) {
        _convergence = conv_fun;
        return solve(function, x_length, x0, x_new);
    }

private:
    IConvergenceChecker *_convergence;
    size_t _max_itr_count;

};

}
