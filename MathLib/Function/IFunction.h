
#pragma once

namespace MathLib
{

template<typename Tval, typename Tpos>
class IFunction
{
public:
	virtual ~IFunction() {};
    virtual Tval eval(const Tpos& x) = 0;
    virtual IFunction<Tval,Tpos>* clone() const = 0;
};

}
