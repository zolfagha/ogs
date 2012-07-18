
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

