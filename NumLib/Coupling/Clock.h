
#pragma once

#include "NumLib/Core/TimeStep.h"
#include "TransientSystems.h"
#include "AsyncPartSolution.h"

namespace NumLib
{

/**
 * \brief Time stepping controller
 */    
class TimeSteppingController
{
public:
    /// set the starting time
	void setBeginning(TimeStep time_begin);
    /// add transient system
	void addTransientSystem(ITransientCoupledProblem &sys);
	/// solve systems until the given time
	void solve(TimeStep time_end);

private:
    TimeStep _time_begin;
    ITransientCoupledProblem* _root_subsystems;
};

}
