
#pragma once

#include "INonlinearSolver.h"

namespace NumLib
{

/**
 * \brief Linear
 */
template <class F_LINEAR>
class Linear : public INonlinearSolver
{
    F_LINEAR* _linear_f;
public:
    explicit Linear(F_LINEAR* linear_f) : _linear_f(linear_f) {};
    virtual ~Linear() {};

    virtual void solve(const VectorType &x_0, VectorType &x_new)
    {
        _linear_f->eval(x_0, x_new);
    }
};

} //end

