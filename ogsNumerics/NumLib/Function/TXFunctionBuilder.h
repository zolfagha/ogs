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

#include "BaseLib/Options.h"
#include "TXFunction.h"

namespace NumLib
{

/**
 * \brief Builder of TX functions
 */
class TXFunctionBuilder
{
public:
    /**
     *
     * @param opDistribution    Option tree for distribution function
     * @return
     */
    static ITXFunction* create(const BaseLib::Options &opDistribution);
};


}
