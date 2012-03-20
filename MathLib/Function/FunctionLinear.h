
#pragma once

#include "IFunction.h"
#include "MathLib/LinearInterpolation.h"

namespace MathLib
{

class FunctionLinear1D : public IFunction<double, double>
{
public:
    FunctionLinear1D(LinearInterpolation* linear) {
        _linear = linear;
    };
    virtual double eval(const double& x)
    {
        return _linear->getValue(x);
    };
    virtual IFunction<double,double>* clone() const
    {
        FunctionLinear1D* obj = new FunctionLinear1D(_linear);
        return obj;
    }
private:
    LinearInterpolation *_linear;
};

} //end
