
#pragma once

#include "NumLib/TimeStepping/TimeStep.h"
#include "NumLib/TimeStepping/ITransientSystem.h"
#include "NamedVariableContainer.h"

namespace NumLib
{

/**
 * \brief Interface class of coupling problems
 */
class ICoupledSystem
{
public:
    virtual int solve() = 0;

    virtual void setParameter(size_t varId, Variable*) = 0;
    virtual Variable* getParameter(size_t varId) const = 0;
    virtual bool check() const= 0;
    virtual size_t getNumberOfInputParameters() const = 0;
    virtual size_t getNumberOfParameters() const = 0;
    virtual size_t getParameterIdForInput(size_t input_id) const = 0;
};

/**
 * \brief Interface class of transient coupling problems
 */
class ITransientCoupledSystem : public ICoupledSystem, public ITransientSystem
{
public:
    int solve(const TimeStep &time) 
    {
        setCurrentTime(time);
        return solve();
    }
private:
    int solve()
    {
        return solveTimeStep(getCurrentTime());
    }
//public:
//    virtual double suggestNext(const TimeStep &time_current) = 0;
//    virtual bool isAwake(const TimeStep &time) = 0;
//    virtual int solveTimeStep(const TimeStep &time) = 0;
//
//    int solve(const TimeStep &time) 
//    {
//        _current_time = time;
//        return solve();
//    }
//
//    const TimeStep& getCurrentTime() const {return _current_time;};
//    void setCurrentTime(const TimeStep &t) {_current_time = t;};
//private:
//    int solve()
//    {
//        return solveTimeStep(_current_time);
//    }
//    TimeStep _current_time;
};

}
