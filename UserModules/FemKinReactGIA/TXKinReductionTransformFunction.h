/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TXKinReductionTransformFunction.h
 *
 * Created on 2012-09-25 by Haibing Shao
 */

#ifndef TX_Kin_Reduction_Transform_Function_H
#define TX_Kin_Reduction_Transform_Function_H

#include "TXPosition.h"
#include "ITXFunction.h"

namespace NumLib
{

/**
 *
 */
template <class T1, class T2, template <typename> class T_OPERATOR>
class TXKinReductionTransformFunction : public ITXFunction
{
public:
    TXKinReductionTransformFunction(T1 &f1, T2 &f2) : _f1(&f1), _f2(&f2) {};
    virtual ~TXKinReductionTransformFunction() {};

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

    virtual TXKinReductionTransformFunction* clone() const
    {
        return new TXKinReductionTransformFunction<T1, T2, T_OPERATOR>(*_f1, *_f2);
    }

private:
    T1 *_f1;
    T2 *_f2;
};

} //end


#endif  // end of ifndef

