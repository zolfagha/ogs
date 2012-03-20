
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

    static double convertToSec(TimeUnit::type type, double v)
    {
        double fac = 1.0;
        switch (type)
        {
        case TimeUnit::Minute:
            fac = 60.0;
        	break;
        case TimeUnit::Hour:
            fac = 3600.0;
        	break;
        case TimeUnit::Day:
            fac = 86400.;
        	break;
        case TimeUnit::Week:
            fac = 86400.*7.0;
        	break;
        case TimeUnit::Year:
            fac = 86400.*365.0;
        	break;
        default:
        	break;
        }
        return v*fac;
    }
};


class ITimeStepFunction
{
public:
    virtual double getBeginning() const = 0;
    virtual double getEnd() const = 0;
    virtual double getPrevious() const = 0;
    virtual double getNext(double t_current) = 0;
    virtual void accept() = 0;
    virtual ITimeStepFunction* clone() = 0;
    virtual ~ITimeStepFunction() {};
};

class AbstractTimeStepFunction : public ITimeStepFunction
{
public:
    AbstractTimeStepFunction(double t0, double tn, TimeUnit::type t_unit)
    {
        _t0 = TimeUnit::convertToSec(t_unit, t0);
        _tn = TimeUnit::convertToSec(t_unit, tn);
        _t_previous = _t_next = _t0;
        _steps = 0;
    };
    virtual ~AbstractTimeStepFunction() {};
    double getBeginning() const {return _t0;};
    double getEnd() const {return _tn;};
    void accept()
    {
        _t_previous = _t_next;
        _steps++;
    }
    double getPrevious() const {return _t_previous;};
    double getNext(double t_current)
    {
        _t_next = suggestNext(t_current);
        return _t_next;
    }

protected:
    virtual double suggestNext(double t_current) = 0;
private:
    double _t0;
    double _tn;
    double _t_previous;
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
        _dt = TimeUnit::convertToSec(t_unit, dt);
    };
    ITimeStepFunction* clone()
    {
        TimeStepFunctionConstant *obj = new TimeStepFunctionConstant(getBeginning(), getEnd(), _dt);
        return obj;
    }

protected:
    double suggestNext(double) 
    {
        double t = getPrevious() + _dt;
        if (t > getEnd())
            return getEnd();
        else
            return t;
    };
};


}
