/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ICoupledProblem.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "NumLib/IOSystem/IIOSystem.h"
#include "NumLib/IOSystem/INamedIO.h"

namespace NumLib
{

/**
 * \brief Interface class of coupling problems
 */
class ICoupledSystem : public IIOSystem, public INamedIO
{
public:
    /// 
    virtual ~ICoupledSystem() {};

    /// solve
    virtual int solve() = 0;

    /// check consistency
    virtual bool check() const= 0;
};

}
