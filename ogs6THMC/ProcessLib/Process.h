/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Process.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "BaseLib/Options.h"
#include "NumLib/TransientCoupling/TransientMonolithicSystem.h"

namespace ProcessLib
{

/**
 * \brief Interface definition of process class
 *
 */
typedef NumLib::AbstractTransientMonolithicSystem Process;

///**
// * \brief Interface definition of process class
// *
// */
//class Process : public NumLib::AbstractTransientMonolithicSystem
//{
//public:
//    /// initialize
//    virtual void initialize(const BaseLib::Options &op) = 0;
//    /// finalize
//    virtual void finalize() = 0;
//};


}

