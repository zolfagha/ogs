
#pragma once

#include "NumLib/TimeStepping/TimeStep.h"

class IOutputTiming
{
public:
    virtual ~IOutputTiming() {};
    virtual bool isActive(const NumLib::TimeStep &current_timestep) = 0;
};

