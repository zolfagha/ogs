
#pragma once

#include <vector>

#include "NumLib/TimeStepping/TimeStep.h"


class IOutputTiming
{
public:
    virtual ~IOutputTiming() {};
    virtual bool isActive(const NumLib::TimeStep &current_timestep) = 0;
};

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

class OutputTimingList : public IOutputTiming
{
public:
    explicit OutputTimingList(std::vector<double> &vec_time)
    : _list_time(vec_time)
    {
    }

    virtual ~OutputTimingList() {};

    virtual bool isActive(const NumLib::TimeStep &current_timestep)
    {
        return (_list_time.end() != std::find(_list_time.begin(), _list_time.end(), current_timestep.getTime()));
    }

private:
    std::vector<double> _list_time;
};

class OutputTimingBuilder
{
public:
    IOutputTiming* create(const std::string &name, size_t n, std::vector<double> *vec_time = 0)
    {
        if (name.compare("STEPS")==0) {
            return new OutputTimingStepPeriodic(n);
        }
        return NULL;
    }
};
