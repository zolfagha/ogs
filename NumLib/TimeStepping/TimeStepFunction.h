
#pragma once

namespace NumLib
{

struct TimeUnit
{
    enum type
    {
        Second,
        Minute,
        Hour,
        Day,
        Week,
        Year
    };
};

double convertTimeUnitToSec(TimeUnit::type type, double v)
{
    double fac = 1.0;
    switch (type)
    {
    case TimeUnit::Minute:
        fac = 60.0;
    case TimeUnit::Hour:
        fac = 3600.0;
    case TimeUnit::Day:
        fac = 86400.;
    case TimeUnit::Week:
        fac = 86400.*7.0;
    case TimeUnit::Year:
        fac = 86400.*365.0;
    }
    return v*fac;
}

class ITimeStepFunction
{
public:
    virtual double getBeginning() const = 0;
    virtual double getEnd() const = 0;
    virtual double next(double t_current) = 0;
    virtual void accept() = 0;
};

class AbstractTimeStepFunction : public ITimeStepFunction
{
public:
    AbstractTimeStepFunction(double t0, double tn, TimeUnit::type t_unit)
    {
        _t0 = convertTimeUnitToSec(t_unit, t0);
        _tn = convertTimeUnitToSec(t_unit, tn);
        _t_current = _t0;
        _t_next = _t_current;
        _steps = 0;
    };
    double getBeginning() const {return _t0;};
    double getEnd() const {return _tn;};
    void accept()
    {
        _t_current = _t_next;
        _steps++;
    }
    double getCurrent() const {return _t_current;};
    double next(double t_current)
    {
        _t_next = suggestNext(t_current);
        return _t_next;
    }
protected:
    virtual double suggestNext(double t_current) = 0;
private:
    double _t0;
    double _tn;
    double _t_current;
    double _t_next;
    size_t _steps;
};

class TimeStepFunctionConstant : public AbstractTimeStepFunction
{
private:
    double _dt;
public:
    TimeStepFunctionConstant(double t0, double tn, double dt, TimeUnit::type t_unit=TimeUnit::Second) : AbstractTimeStepFunction(t0, tn, t_unit)
    {
        _dt = convertTimeUnitToSec(t_unit, dt);
    };

protected:
    double suggestNext(double) 
    {
        double t = getCurrent() + _dt;
        if (t > getEnd())
            return getEnd();
        else
            return t;
    };
};


}