
#pragma once

#include "BaseLib/Options.h"
#include "NumLib/TransientCoupling/TransientMonolithicSystem.h"

namespace ProcessLib
{

//typedef NumLib::AbstractTransientMonolithicSystem Process;

class Process : public NumLib::AbstractTransientMonolithicSystem
{
public:
	virtual void initialize(const BaseLib::Options &op) = 0;
	virtual void finalize() = 0;
};


}

