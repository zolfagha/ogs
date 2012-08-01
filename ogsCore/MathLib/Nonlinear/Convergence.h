
#pragma once

#include <cmath>
#include <algorithm>

#include "logog/include/logog.hpp"

#include "MathLib/LinAlg/VectorNorms.h"


namespace MathLib
{

/**
 * \brief Convergence check for Newton iterations
 */
template<class T_D0, class T_ERROR>
class NRCheckConvergence
{
public:
    NRCheckConvergence() : _tolerance(1.e-6), _error(.0) {};
    explicit NRCheckConvergence(double err) : _tolerance(err), _error(.0) {};

    bool check(T_D0* r, T_D0* dx, T_D0* x_new)
    {
        _error = calc.error(r, dx, x_new);
        return (fabs(_error) < _tolerance);
    }

    double getError() const {return _error;};

    double getTolerance() const {return _tolerance;};

private:
    double _tolerance;
    double _error;
    T_ERROR calc;
};


/**
 * \brief Check norm1 of dx
 */
class NRErrorNorm1DX
{
public:
    template<class T_D0>
    inline double error(T_D0*, T_D0* dx, T_D0*)
    {
        if (dx==0) return .0;
        return norm1(*dx, dx->size());
    }
};

template<>
inline double NRErrorNorm1DX::error(double*, double* dx, double*)
{
    return fabs(*dx);
}

/**
 * \brief Check norm1 of residual
 */
class NRErrorNorm1Residual
{
public:
    template<class T_D0>
    inline double error(T_D0* r, T_D0*, T_D0*)
    {
        if (r==0) return .0;
        return norm1(*r, r->size());
    }
};

template<>
inline double NRErrorNorm1Residual::error(double* r, double*, double*)
{
    return fabs(*r);
}

/**
 * \brief Check norm-max of residual and dx. Relative error is used for dx.
 */
class NRErrorAbsResMNormOrRelDxMNorm
{
private:
    size_t _itr_count;
public:
    NRErrorAbsResMNormOrRelDxMNorm() : _itr_count(0) {};

    template<class T_D0>
    inline double error(T_D0* r, T_D0* dx, T_D0* x)
    {
        double abs_mnorm_r = norm_max(*r, r->size());
        double abs_mnorm_x = norm_max(*x, x->size());
        double abs_mnrom_dx = norm_max(*dx, dx->size());
        return calc(abs_mnorm_r, abs_mnorm_x, abs_mnrom_dx);
    }

private:
    inline double calc(double abs_mnorm_r, double abs_mnorm_x, double abs_mnrom_dx)
    {
        double rel_mnorm_dx = .0;
        if (abs_mnorm_x!=.0) rel_mnorm_dx = abs_mnrom_dx/abs_mnorm_x;

        INFO("-> %d: ||r||_inf=%1.3e, ||dx||_inf=%1.3e, ||dx||/||x||=%1.3e", _itr_count, abs_mnorm_r, abs_mnrom_dx, rel_mnorm_dx);

        _itr_count++;

        if (abs_mnorm_x == .0) {
            return abs_mnorm_r;
        } else {
            return std::max(abs_mnorm_r, rel_mnorm_dx);
        }
    }
};

template<>
inline double NRErrorAbsResMNormOrRelDxMNorm::error(double* r, double* dx, double* x)
{
    double abs_mnorm_r = fabs(*r);
    double abs_mnorm_x = fabs(*x);
    double abs_mnrom_dx = fabs(*dx);
    return calc(abs_mnorm_r, abs_mnorm_x, abs_mnrom_dx);
}

}
