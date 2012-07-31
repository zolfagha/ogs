
#pragma once

#include "BaseLib/CodingTools.h"
#include "NumLib/TimeStepping/TimeSteppingController.h"
#include "OutputController.h"

namespace ogs6
{

class TimeSteppingControllerWithOutput : public NumLib::TimeSteppingController
{
public:
    TimeSteppingControllerWithOutput(OutputController* p) : _out_controller(p) {};

    virtual ~TimeSteppingControllerWithOutput()
    {
    };

protected:
    virtual void doSomethingAfterTimeStepAccepted(const NumLib::TimeStep &t_current) const;

private:
    OutputController* _out_controller;
};

} //ogs6