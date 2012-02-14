
#pragma once

#include "TimeStep.h"
#include "ITransientProblem.h"

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
	void addTransientSystem(ITransientProblem &sys);
	/// solve systems until the given time
	void solve(double time_end);

private:
    double _time_begin;
    ITransientProblem* _root_subsystems;
};

}
