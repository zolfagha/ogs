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

class TXFunctionBuilder
{
public:
    ITXFunction* create(const std::string &name, double v)
    {
        if (name.compare("CONSTANT")==0) {
            return new TXFunctionConstant(v);
        }
        return NULL;
    }
};


}
