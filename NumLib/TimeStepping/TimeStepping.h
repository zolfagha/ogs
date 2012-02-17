
#pragma once

#include <vector>

#include "TimeStep.h"

namespace NumLib
{
class ITimeStepping
{
public:
    virtual double next() const = 0;
    virtual double begin() const = 0;
    virtual double end() const = 0;
};

class FixedTimeStepping : public ITimeStepping
{
private:
    std::vector<double> _list_timesteps;
    size_t _count;
public:
    FixedTimeStepping() : _count(0) {};
    void setTimeSteps(std::vector<double> &v)
    {
        _list_timesteps.assign(v.begin(), v.end());
    }
    double begin() const 
    {
        return *_list_timesteps.begin();
    }
    double end() const
    {
        return _list_timesteps.back();
    }
    double next() const {
        return _list_timesteps[_count];
    }

    size_t getNumberOfTimeStep() const {return _list_timesteps.size();};


    void accepted() { _count++;};
};

}
