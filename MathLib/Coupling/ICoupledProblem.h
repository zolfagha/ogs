
#pragma once

#include "MathLib/Parameter/IIOSystem.h"
#include "MathLib/Parameter/NamedIO.h"

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
