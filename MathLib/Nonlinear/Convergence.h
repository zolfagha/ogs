
#pragma once

#include <cmath>

namespace MathLib {

class VectorNormCalculater
{
public:
    static double calculate(double *x, size_t length, int norm_type);
};

class IConvergenceChecker
{
public:
    virtual bool isConverged(double *x) = 0;
};

class ConvergenceChecker : public IConvergenceChecker
{
public:
    ConvergenceChecker() {
        _normType = 1;
        _tolerance = 1.e-10;
    };

    ~ConvergenceChecker() {};

    virtual bool isConverged(double *x, size_t length) {
        double norm = VectorNormCalculater::calculate(x, length, _normType);
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
