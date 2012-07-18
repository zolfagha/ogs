
#pragma once

#include "NumLib/TimeStepping/TimeStep.h"
#include "IOutputTiming.h"

class OutputTimingStepPeriodic : public IOutputTiming
{
public:
    explicit OutputTimingStepPeriodic(size_t step)
    : _steps(step)
    {
    }

    virtual ~OutputTimingStepPeriodic() {};

    virtual bool isActive(const NumLib::TimeStep &current_timestep)
    {
        return (current_timestep.getTimeStepCount() % _steps == 0);
    }

private:
    size_t _steps;
};
