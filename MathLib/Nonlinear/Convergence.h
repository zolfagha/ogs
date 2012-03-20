
#pragma once

#include <cmath>

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

}
