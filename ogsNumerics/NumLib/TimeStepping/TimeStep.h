
#pragma once

#include "BaseLib/CodingTools.h"

namespace NumLib
{

/**
 * \brief Time step class
 */    
class TimeStep
{
public:
    TimeStep() : _time_current(.0), _dt(.0), _time_stepping_count(0) {};
    TimeStep(double t) : _time_current(t), _dt(.0), _time_stepping_count(0) {};
    TimeStep(double t, double dt) : _time_current(t), _dt(dt), _time_stepping_count(0) {};
    TimeStep(const TimeStep &prev, double t) 
    {
        _time_current = t;
        _dt = t - prev._time_current;
        _time_stepping_count = prev._time_stepping_count+1;
    };

    bool operator<(const TimeStep &t) const
    {
        return (_time_current < t._time_current);
    }

    size_t getTimeStepCount() const {return _time_stepping_count;};
    double getTime() const {return _time_current;};
    double getTimeStepSize() const {return _dt;};

    void setTime(double t) {_time_current = t;};
    void setTimeStepSize(double dt) {_dt = dt;};
    void accept(double t) 
    {
        _time_current = t;
        _dt = .0;
        _time_stepping_count++;
    };

    void assign(const TimeStep& src)
    {
        _time_current = src._time_current;
        _dt = src._dt;
        _time_stepping_count = src._time_stepping_count;
    }

private:
    double _time_current;
    double _dt;
    size_t _time_stepping_count;

    DISALLOW_COPY_AND_ASSIGN(TimeStep);
    //TimeStep(const TimeStep &prev) 
    //{
    //    _time_current = prev._time_current;
    //    _dt = prev._dt;
    //    _time_stepping_count = prev._time_stepping_count;
    //};
};

}