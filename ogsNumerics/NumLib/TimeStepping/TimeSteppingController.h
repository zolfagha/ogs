
#pragma once

#include "TimeStep.h"
#include "ITransientSystem.h"

namespace NumLib
{


/**
 * \brief Time stepping controller
 */    
class TimeSteppingController
{
public:
    /// set the starting time
	void setBeginning(double time_begin);
    /// add transient system
	void addTransientSystem(ITransientSystem &sys);
	/// solve systems until the given time
	size_t solve(double time_end);

private:
    double _time_begin;
    ITransientSystem* _root_subsystems;
};

}
