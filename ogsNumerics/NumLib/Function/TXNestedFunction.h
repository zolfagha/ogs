/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 */

#pragma once

#include "TXPosition.h"
#include "ITXFunction.h"

namespace NumLib
{

/**
 * \brief Temporal-spatial function made of two functions
 *
 */
template <class T1, class T2>
class TXNestedFunction : public ITXFunction
{
public:
    /**
     *
     * @param f1
     * @param f2
     */
	TXNestedFunction(T1** f1, T2* f2) : _f1(f1), _f2(f2)
    {
        isConst((*_f1)->isConst());
        isTemporallyConst((*_f1)->isTemporallyConst());
        isSpatiallyConst((*_f1)->isSpatiallyConst());
    };

    /**
     * Copy constructor
     * @param src
     */
	TXNestedFunction(const TXNestedFunction &src)
    : _f1(src._f1), _f2(src._f2->clone())
    {
        isConst((*_f1)->isConst());
        isTemporallyConst((*_f1)->isTemporallyConst());
        isSpatiallyConst((*_f1)->isSpatiallyConst());
    }

    /**
     *
     */
    virtual ~TXNestedFunction()
    {
        BaseLib::releaseObject(_f2);
    };

    /**
     *
     * @param x
     * @param val
     */
    virtual void eval(const TXPosition x, DataType &val) const OGS_DECL_OVERRIDE
    {
        DataType v1;
        (*_f1)->eval(x, v1);
        _f2->eval(v1, val);
    }

    /**
     *
     * @param x
     * @param v
     */
    virtual void eval(const TXPosition x, double &val) const OGS_DECL_OVERRIDE
    {
        double v1;
        (*_f1)->eval(x, v1);
        _f2->eval(v1, val);
    }

    /**
     *
     * @param x
     * @param val
     */
    virtual void eval(const double* x, double &val) const OGS_DECL_OVERRIDE
    {
        double v1;
        (*_f1)->eval(x, v1);
        _f2->eval(v1, val);
    }

    /**
     *
     * @return
     */
    virtual TXNestedFunction* clone() const OGS_DECL_OVERRIDE
    {
        return new TXNestedFunction<T1, T2>(*this);
    }

private:
    T1** _f1;
    T2* _f2;
};

} //end


