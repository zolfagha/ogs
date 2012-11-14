/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TXFunctionAnalytical.h
 *
 * Created on 2012-11-14 by Norihiro Watanabe
 */

#pragma once

#include <string>
#include <vector>
#include "TXPosition.h"
#include "ITXFunction.h"

namespace NumLib
{

/**
 * \brief Constant value
 */
class TXFunctionAnalytical : public ITXFunction
{
public:
    TXFunctionAnalytical(const std::string &str_expression)
    {
        ITXFunction::isConst(true);
        ITXFunction::isTemporallyConst(true);
        ITXFunction::isSpatiallyConst(false);
    };

    virtual ~TXFunctionAnalytical() {};

    virtual void eval(const TXPosition x, double &val) const
    {
    }

    virtual TXFunctionAnalytical* clone() const
    {
        return new TXFunctionAnalytical(_exp);
    }

private:
    std::string _exp;
};

} //end


