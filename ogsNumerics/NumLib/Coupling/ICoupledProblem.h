
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
