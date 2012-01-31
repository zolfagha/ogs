
#pragma once

#include "NumLib/Core/TimeStep.h"

namespace NumLib
{

class ITransientSystem
{
public:
	virtual TimeStep suggestNext(TimeStep time_current) = 0;
	virtual bool solveNextStep(TimeStep time) = 0;
	virtual bool isAwake(TimeStep time) = 0;
};



}
