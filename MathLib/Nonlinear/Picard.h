
#pragma once

#include "Nonlinear/Convergence.h"

namespace MathLib {

class Picard
{
public:
    Picard() {
        _max_itr_count = 100;
        _convergence = new ConvergenceChecker();
    }
    virtual ~Picard() {
        if (_convergence)
            delete _convergence;
    }

    //TODO should user prepare memory for x_new?
    int solve(Function *function, size_t x_length, double* x0, double *x_new) {

        double *x_prev = x0;
        for (size_t i=0; i<_max_itr_count; i++) {
            x_new = function->update(x_prev);
            if (_convergence->isConverged(x_new)) 
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
