
#pragma once

#include "MathLib/IOSystem/IIOSystem.h"
#include "MathLib/IOSystem/INamedIO.h"

namespace MathLib
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
