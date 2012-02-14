
#pragma once

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

    size_t getTimeStep() const {return _time_stepping_count;};
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

private:
    double _time_current;
    double _dt;
    size_t _time_stepping_count;

    TimeStep(const TimeStep &prev) 
    {
        _time_current = prev._time_current;
        _dt = prev._dt;
        _time_stepping_count = prev._time_stepping_count;
    };
};

}