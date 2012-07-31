
#include "TimeSteppingControllerWithOutput.h"

namespace ogs6
{

void TimeSteppingControllerWithOutput::doSomethingAfterTimeStepAccepted(const NumLib::TimeStep &t_current) const
{
    _out_controller->outputData(t_current);
}

}
