
#pragma once

#include "NumLib/Core/TimeStep.h"
#include "NamedVariableContainer.h"

namespace NumLib
{

/**
 * \brief Interface class of coupling problems
 */
class ICoupledProblem
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
class ITransientCoupledProblem : public ICoupledProblem
{
public:
    virtual TimeStep suggestNext(TimeStep time_current) = 0;
    virtual bool isAwake(TimeStep time) = 0;
    virtual int solveTimeStep(TimeStep time) = 0;

    int solve(TimeStep time) {
        _current_time = time;
        return solve();
    }



    TimeStep getCurrentTime() const {return _current_time;};
    void setCurrentTime(TimeStep t) {_current_time = t;};
private:
    int solve()
    {
        return solveTimeStep(_current_time);
    }
    TimeStep _current_time;
};

}
