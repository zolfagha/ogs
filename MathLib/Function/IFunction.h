
#pragma once

#include "MathLib/Vector.h"
#include "MathLib/LinAlg/Dense/Matrix.h"
#include "MathLib/IClonable.h"

namespace MathLib
{

class IFunction : public IClonable
{
public:
	virtual ~IFunction() {};
};

template<typename Tpos, typename Tval>
class TemplateFunction : public IFunction
{
public:
	virtual ~TemplateFunction() {};
    virtual void eval(const Tpos &x, Tval &val) = 0;
//    virtual TemplateFunction<Tpos,Tval>* clone() const = 0;
};

template<typename Tpos, typename Tval, class T_FUNCTION>
class AbstractDefaultCloneFunction : public TemplateFunction<Tpos, Tval>
{
public:
	virtual ~AbstractDefaultCloneFunction() {};
    virtual IFunction* clone() const
	{
    	return new T_FUNCTION();
	}
};

typedef double* SpatialPosition;

template<typename Tval>
class SpatialFunction : public TemplateFunction<SpatialPosition, Tval>
{
public:
	virtual ~SpatialFunction() {};
};

typedef SpatialFunction<double> SpatialFunctionScalar;
typedef SpatialFunction<MathLib::Vector> SpatialFunctionVector;
typedef SpatialFunction<MathLib::Matrix<double> > SpatialFunctionTensor;

template <typename Tval>
class Multiplication
{
public:
    void doit(const Tval &v1, const Tval &v2, Tval &val)
    {
        val = v1 * v2;
    }
};

template <typename Tval, class T1, class T2, template <typename> class T_OPERATOR>
class SpatialCompositFunction : public SpatialFunction<Tval>
{
    T1 *_f1;
    T2 *_f2;
public:
    SpatialCompositFunction(T1 &f1, T2 &f2) : _f1(&f1), _f2(&f2) {};
    virtual ~SpatialCompositFunction() {};
    virtual void eval(const SpatialPosition &x, Tval &val)
    {
        Tval v1, v2;
        _f1->eval(x, v1);
        _f2->eval(x, v2);
        T_OPERATOR<Tval> op;
        op.doit(v1, v2, val);
    }
    virtual SpatialCompositFunction* clone() const
    {
        return new SpatialCompositFunction<Tval, T1, T2, T_OPERATOR>(*_f1, *_f2);
    }


};

}
