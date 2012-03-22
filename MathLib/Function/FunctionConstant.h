
#pragma once

#include "IFunction.h"

namespace MathLib
{

template<typename Tpos, typename Tval>
class FunctionConstant : public IFunction<Tpos,Tval>
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

    virtual IFunction<Tpos,Tval>* clone() const
    {
        IFunction<Tpos,Tval>* obj = new FunctionConstant(_v);
        return obj;
    }
private:
    Tval _v;
};

} //end
