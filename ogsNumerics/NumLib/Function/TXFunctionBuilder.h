/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TXFunctionBuilder.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include <string>
#include "TXFunction.h"

namespace NumLib
{

/**
 * \brief Builder of TX functions
 */
class TXFunctionBuilder
{
public:
    ITXFunction* create(const TXFunctionType::type f_type, double v)
    {
        switch (f_type) {
        case TXFunctionType::CONSTANT:
            return new TXFunctionConstant(v);
        case TXFunctionType::LINEAR:
        default:
            break;
        }

        return NULL;
    }

    ITXFunction* create(const std::string &name, double v)
    {
        return create(convertStringToTXFunctionType(name), v);
    }
};


}
