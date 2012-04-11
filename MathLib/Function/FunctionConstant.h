
#pragma once

#include "IFunction.h"

namespace MathLib
{

template<typename Tpos, typename Tval>
class FunctionConstant : public TemplateFunction<Tpos,Tval>
{
public:
    FunctionConstant(const Tval &v)
    {
        _v = v;
    };

    virtual ~FunctionConstant() {};

    virtual void eval(const Tpos&, Tval &v)
    {
        v = _v;
    };

    virtual void eval(const Tpos&, Tval &v) const
    {
        v = _v;
    };

    virtual void eval(Tval &v) const
    {
        v = _v;
    };

    virtual TemplateFunction<Tpos,Tval>* clone() const
    {
    	return new FunctionConstant(_v);
    }
private:
    Tval _v;
};

template<typename Tval>
class SpatialFunctionConstant : public SpatialFunction<Tval>
{
public:
    SpatialFunctionConstant(const Tval &v)
    {
        _v = v;
    };

    virtual ~SpatialFunctionConstant() {};

    virtual void eval(const SpatialPosition&, Tval &v)
    {
        v = _v;
    };

    virtual void eval(const SpatialPosition&, Tval &v) const
    {
        v = _v;
    };

    virtual void eval(Tval &v) const
    {
        v = _v;
    };

    virtual SpatialFunctionConstant<Tval>* clone() const
    {
        return new SpatialFunctionConstant<Tval>(_v);
    }
private:
    Tval _v;
};
} //end
