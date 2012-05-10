
#pragma once

#include "NumLib/Function/IFunction.h"
#include "GeoLib/Core/Point.h"

namespace NumLib
{

struct TXPosition
{
	TXPosition(double t_, double* x_) : t(t_), x(x_) {};
	TXPosition(double* x_) : x(x_) {};
	TXPosition(const GeoLib::Point* p) : x(p->getData()) {};

	double t;
	const double* x;
};

class ITXValue
{
public:
	virtual ~ITXValue() {};
	virtual bool isScalar() const = 0;
	virtual bool isVector() const = 0;
	virtual bool isTensor() const = 0;
};

class TXValueScalar
{
public:
	virtual ~TXValueScalar() {};
	bool isScalar() const {return true;};
	bool isVector() const {return false;};
	bool isTensor() const {return false;};
};

class ITXFunction
{
public:
	virtual ~ITXFunction() {};
	virtual void eval(const TXPosition x, ITXValue &v) = 0;
	virtual bool isConst() = 0;
	virtual bool isTemporallyConst() = 0;
	virtual bool isSpatiallyConst() = 0;
};

class ITXFunctionScalar : public ITXFunction
{
public:
	virtual ~ITXFunctionScalar() {};
//	void eval(const TXPosition x, ITXValue &v)
//	{
//		eval(x, (TXValueScalar)v);
//	}
//	virtual void eval(const TXPosition x, TXValueScalar &v) = 0;
	virtual void eval(const TXPosition x, double &v) = 0;
	virtual bool isConst() = 0;
	virtual bool isTemporallyConst() = 0;
	virtual bool isSpatiallyConst() = 0;
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

typedef TXFunctionConstant<double, ITXFunctionScalar> TXFunctionConstantScalar;

} //end
