
#pragma once

#include "NumLib/Core/TimeStep.h"
#include "TransientSystems.h"
#include "AsyncPartSolution.h"

namespace NumLib
{

class Clock
{
protected:
	TimeStep _time_begin;
	AsyncPartSolution rootAsyncPartSolution;

public:
	void setBeginning(TimeStep time);

	void addTransientSystem(ITransientSystem *system);
	
	void moveForwardUntill(TimeStep time_end);
};

}
