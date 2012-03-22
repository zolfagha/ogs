
#pragma once

namespace MathLib
{

template<typename Tpos, typename Tval>
class IFunction
{
public:
	virtual ~IFunction() {};
    virtual void eval(const Tpos &x, Tval &val) = 0;
    virtual IFunction<Tpos,Tval>* clone() const = 0;
};

template<typename Tpos, typename Tval, class T_FUNCTION>
class AbstractDefaultCloneFunction : public IFunction<Tpos, Tval>
{
public:
	virtual ~AbstractDefaultCloneFunction() {};
    virtual IFunction<Tpos,Tval>* clone() const
	{
    	return new T_FUNCTION();
	}
};

}
