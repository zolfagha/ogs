
#pragma once

#include <vector>

#include "NumLib/TimeStepping/ITransientSystem.h"
#include "NumLib/TimeStepping/TimeStepFunction.h"

#include "ISolutionAlgorithm.h"

namespace SolutionLib
{

/**
 * \brief Abstract class for time-stepping method
 */
class AbstractTimeSteppingAlgorithm : public ISolutionAlgorithm, public NumLib::ITransientSystem
{
public:
    /// @param tim Time step function
    AbstractTimeSteppingAlgorithm(NumLib::ITimeStepFunction &tim) : _tim(&tim) {};

    /// destructor
    virtual ~AbstractTimeSteppingAlgorithm() {};

    /// get the time step function
    /// @return Time step function
    NumLib::ITimeStepFunction* getTimeStepFunction() const {return _tim;};

    /// suggest the next time step
    /// @param time_currrent current time step object
    /// @return the next time
    double suggestNext(const NumLib::TimeStep &time_current)
    {
        return _tim->getNext(time_current.getTime());
    }

    /// return if the next time step is the same as the given time
    /// @param time the time step object
    /// @bool true if this process is active with the given time
    bool isAwake(const NumLib::TimeStep &time)
    {
        return (time.getTime()==_tim->getNext(time.getTime()));
    }

    virtual void accept(const NumLib::TimeStep &)
    {
        _tim->accept();
    };


private:
    NumLib::ITimeStepFunction* _tim;
};


}
