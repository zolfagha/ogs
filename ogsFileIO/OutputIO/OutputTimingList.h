
#pragma once

#include <vector>

#include "NumLib/TimeStepping/TimeStep.h"
#include "IOutputTiming.h"

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
