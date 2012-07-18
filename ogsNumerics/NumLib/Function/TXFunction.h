
#pragma once

#include "NumLib/DataType.h"
#include "NumLib/Function/IFunction.h"
#include "TXPosition.h"

namespace NumLib
{

/**
 * \brief Interface of any functions in space-time domain
 *
 * This class aims to be an abstract of spatially and temporally distributed data such as
 * physical quantity (e.g. head, stress) and material property (e.g. permeability).
 *
 * TXFunction
 * - is evaluated at particular position in space-time domain
 * - returns scalar or vector value
 * - has some attributes (e.g constant)
 *
 */
class ITXFunction : public IClonable
{
public:
    typedef LocalMatrix DataType;
    virtual ~ITXFunction() {};

    /// evaluate this function at the given position and return vector data

    /// \param x  position in space and time
    /// \param v  evaluated vector
    virtual void eval(const TXPosition /*x*/, DataType &/*v*/) const {};

    /// evaluate this function at the given position and return scalar value.
    ///
    /// Default behavior of this function is to return the 1st component of the vector.
    /// \param x  position in space and time
    /// \param v  evaluated scalar
    virtual void eval(const TXPosition x, double &v) const
    {
        DataType tmp(1,1);
        this->eval(x, tmp);
        v = tmp(0);
    }

    virtual void eval(const double* x, DataType &v) const
    {
        TXPosition pos(x);
        this->eval(pos, v);
    }

    virtual void eval(const double* x, double &v) const
    {
        TXPosition pos(x);
        this->eval(pos, v);
    }

    virtual ITXFunction* clone() const = 0;

    ///
    virtual bool isConst() const {return false;};
    ///
    virtual bool isTemporallyConst() const {return false;};
    ///
    virtual bool isSpatiallyConst() const {return false;};
};

/**
 * \brief Constant value
 */
class TXFunctionConstant : public ITXFunction
{
public:
    TXFunctionConstant(double val) : _vec(1,1) { _vec(0,0) = val;};
    TXFunctionConstant(const DataType &val) : _vec(val) {};

    virtual ~TXFunctionConstant() {};

    virtual void eval(const TXPosition /*x*/, DataType &v) const { v = _vec;}
    void eval(double &v) const { v = _vec(0,0); };
    void eval(DataType &v) const { v = _vec; };

    virtual TXFunctionConstant* clone() const { return new TXFunctionConstant(_vec); }

    virtual bool isConst() const {return true;};
    virtual bool isTemporallyConst() const {return true;};
    virtual bool isSpatiallyConst() const {return true;};

private:
    DataType _vec;
};


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
template <class T1, class T2, template <typename> class T_OPERATOR>
class TXCompositFunction : public ITXFunction
{
public:
    TXCompositFunction(T1 &f1, T2 &f2) : _f1(&f1), _f2(&f2) {};
    virtual ~TXCompositFunction() {};

    virtual void eval(const TXPosition &x, DataType &val) const
    {
        DataType v1, v2;
        _f1->eval(x, v1);
        _f2->eval(x, v2);
        T_OPERATOR<DataType> op;
        op.doit(v1, v2, val);
    }

    virtual void eval(const double* x, double &val) const
    {
        double v1, v2;
        _f1->eval(x, v1);
        _f2->eval(x, v2);
        T_OPERATOR<double> op;
        op.doit(v1, v2, val);
    }

    virtual TXCompositFunction* clone() const
    {
        return new TXCompositFunction<T1, T2, T_OPERATOR>(*_f1, *_f2);
    }

private:
    T1 *_f1;
    T2 *_f2;
};
} //end


