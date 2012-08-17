/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file MyConvergenceCheckerFactory.h
 *
 * Created on 2012-07-17 by Norihiro Watanabe
 */

#pragma once

#include <string>

#include "NumLib/Function/DiscreteDataConvergenceCheck.h"

class MyConvergenceCheckerFactory
{
public:
    NumLib::IConvergenceCheck* create(const std::string &name)
    {
        if (name.compare("FemFunctionConvergenceCheck")==0) {
            return new NumLib::DiscreteDataConvergenceCheck();
        }
        return NULL;
    };
};
