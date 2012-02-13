
#pragma once

#include "MathLib/LinearInterpolation.h"

namespace MathLib
{

template<typename Tval, typename Tpos>
class IFunction
{
public:
    virtual Tval eval(const Tpos& x) = 0;
    virtual IFunction<Tval,Tpos>* clone() const = 0;
};

template<typename Tval, typename Tpos>
class FunctionConstant : public IFunction<Tval, Tpos>
{
public:
    FunctionConstant(const Tval &v) 
    {
        _v = v;
    };

    virtual Tval eval(const Tpos& x) 
    {
        return _v;
    };

    virtual IFunction<Tval,Tpos>* clone() const
    {
        IFunction<Tval,Tpos>* obj = new FunctionConstant(_v);
        return obj;
    }
private:
    Tval _v;
};

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

}
