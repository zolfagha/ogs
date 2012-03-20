
#pragma once

#include "IFunction.h"

namespace MathLib
{

template<typename Tval, typename Tpos>
class FunctionConstant : public IFunction<Tval, Tpos>
{
public:
    FunctionConstant(const Tval &v)
    {
        _v = v;
    };

    virtual ~FunctionConstant() {};

    virtual Tval eval(const Tpos&)
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

} //end
