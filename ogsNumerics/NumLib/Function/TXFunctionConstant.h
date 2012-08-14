/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TXFunctionConstant.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "TXPosition.h"
#include "ITXFunction.h"

namespace NumLib
{

/**
 * \brief Constant value
 */
class TXFunctionConstant : public ITXFunction
{
public:
    explicit TXFunctionConstant(double val) : _vec(1,1)
    {
        _vec(0,0) = val;
        ITXFunction::isConst(true);
        ITXFunction::isTemporallyConst(true);
        ITXFunction::isSpatiallyConst(true);
    };
    explicit TXFunctionConstant(const DataType &val) : _vec(val) {};

    virtual ~TXFunctionConstant() {};

    virtual void eval(const TXPosition /*x*/, DataType &v) const { v = _vec;}
    void eval(double &v) const { v = _vec(0,0); };
    void eval(DataType &v) const { v = _vec; };

    virtual TXFunctionConstant* clone() const { return new TXFunctionConstant(_vec); }

    //virtual bool isConst() const {return true;};
    //virtual bool isTemporallyConst() const {return true;};
    //virtual bool isSpatiallyConst() const {return true;};

private:
    DataType _vec;
};

} //end


