
#pragma once

#include "BaseLib/Options.h"
#include "TimeStep.h"

namespace NumLib
{

/**
 * \brief Interface class of transient coupling problems
 */
class ITransientSystem
{
public:
    virtual ~ITransientSystem() {};

    /// initialize
    virtual void initialize(const BaseLib::Options &op) {};
    /// finalize
    virtual void finalize() {};

    virtual double suggestNext(const TimeStep &time_current) = 0;
    virtual bool isAwake(const TimeStep &time) = 0;
    virtual int solveTimeStep(const TimeStep &time) = 0;
    const TimeStep& getCurrentTime() const {return *_current_time;};
    void setCurrentTime(const TimeStep &t) {_current_time = const_cast<TimeStep*>(&t);};
    virtual void accept(const TimeStep &time) = 0;

private:
    TimeStep* _current_time;
};




}
