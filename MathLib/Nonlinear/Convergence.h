
#pragma once

#include <cmath>
#include "MathLib/LinAlg/VectorNorms.h"


namespace MathLib
{

class VectorNormCalculater
{
public:
    static double calculate(double *x, size_t length, int norm_type);
};

class IConvergenceChecker
{
public:
	virtual ~IConvergenceChecker() {};
    virtual bool isConverged(double *x, size_t length) = 0;
};

class ConvergenceChecker // : public IConvergenceChecker
{
public:
    ConvergenceChecker() {
        _normType = 1;
        _tolerance = 1.e-10;
    };

    virtual ~ConvergenceChecker() {};

    template<class T_VECTOR>
    bool isConverged(T_VECTOR &x) {
        //double norm = VectorNormCalculater::calculate(x, length, _normType);
    	double norm = 1;
        bool converged = (fabs(norm)<_tolerance);
        return converged;
    }
    void initialize(int norm_type, double tolerance) {
        _normType = norm_type;
        _tolerance = tolerance;
    }

private:
    int _normType;
    double _tolerance;
};

/**
 * \brief Convergence check for Newton iterations
 */
template<class T_D0, class T_ERROR>
class NRCheckConvergence
{
private:
	double _err;
public:
	NRCheckConvergence() : _err(1.e-6) {};
	NRCheckConvergence(double err) : _err(err) {};
	bool check(T_D0* r, T_D0* dx, T_D0* x_new)
	{
		T_ERROR calc;
		double error = calc.error(r, dx, x_new);
		return (fabs(error) < _err);
	}
};

/**
 * \brief Use norm1 of DX as iteration errors
 */
class NRErrorNorm1DX
{
public:
	template<class T_D0>
	inline double error(T_D0* r, T_D0* dx, T_D0* x_new)
	{
		if (dx==0) return .0;
		return norm1(*dx, dx->size());
	}
};

template<>
inline double NRErrorNorm1DX::error(double* r, double* dx, double* x_new)
{
	return *dx;
}

}
