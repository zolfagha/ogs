
#pragma once

#include "NumLib/Function/IFunction.h"
#include "GeoLib/Core/Point.h"

namespace NumLib
{

/**
 * \brief Position in space-time domain
 */
class TXPosition
{
public:
	TXPosition(double t_, double* x_) : _t(t_), _x(x_) {};
	TXPosition(double* x_) : _x(x_) {};
	TXPosition(const GeoLib::Point* p) : _x(p->getData()) {};

	double getTime() const {return _t;};
	const double* getSpace() const {return _x;};
private:
	double _t;
	const double* _x;
};

/**
 * \brief Value
 */
class TXValue
{
public:


private:
};

/**
 * \brief Interface to functions in space-time domain
 */
class ITXFunction : public IClonable
{
public:
	virtual ~ITXFunction() {};
	virtual void eval(const TXPosition x, TXValue &v) = 0;
	virtual void eval(const TXPosition x, double &v) const = 0;
	virtual bool isConst() = 0;
	virtual bool isTemporallyConst() = 0;
	virtual bool isSpatiallyConst() = 0;
};



class ITXValue
{
public:
	virtual ~ITXValue() {};
	virtual bool isScalar() const = 0;
	virtual bool isVector() const = 0;
	virtual bool isTensor() const = 0;
};


template<typename Tval, typename T_BASE>
class TXFunctionConstant : public T_BASE
{
public:
	TXFunctionConstant(const Tval v)
    {
        _v = v;
    };

    virtual ~TXFunctionConstant() {};

	virtual void eval(const TXPosition x, ITXValue &v)
	{
	}

	void eval(const TXPosition x, Tval &v)
	{
		v = _v;
	}

    void eval(Tval &v) const
    {
        v = _v;
    };

	bool isConst() {return true;};
	bool isTemporallyConst() {return true;};
	bool isSpatiallyConst() {return true;};

private:
    Tval _v;
};

typedef TXFunctionConstant<double, ITXFunction> TXFunctionConstantScalar;

/**
 *
 */
template <typename Tval>
class Multiplication
{
public:
    void doit(const Tval &v1, const Tval &v2, Tval &val)
    {
        val = v1 * v2;
    }
};

/**
 *
 */
template <typename Tval, class T1, class T2, template <typename> class T_OPERATOR>
class TXCompositFunction : public ITXFunction
{
    T1 *_f1;
    T2 *_f2;
public:
    TXCompositFunction(T1 &f1, T2 &f2) : _f1(&f1), _f2(&f2) {};
    virtual ~TXCompositFunction() {};
    virtual void eval(const TXPosition &x, Tval &val)
    {
        Tval v1, v2;
        _f1->eval(x, v1);
        _f2->eval(x, v2);
        T_OPERATOR<Tval> op;
        op.doit(v1, v2, val);
    }
    virtual TXCompositFunction* clone() const
    {
        return new TXCompositFunction<Tval, T1, T2, T_OPERATOR>(*_f1, *_f2);
    }


};
} //end
