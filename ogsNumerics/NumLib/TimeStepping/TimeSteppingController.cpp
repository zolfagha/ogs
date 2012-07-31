
#include "TimeSteppingController.h"

#include <iostream>

#include "logog/include/logog.hpp"

namespace NumLib
{

size_t TimeSteppingController::solve(double time_end) 
{
    TimeStep time_current(_time_begin);
    INFO("\n\n***Start time steps with t0=%d s\n", _time_begin);

    while (time_current.getTime()<time_end) {
        double time_next = _root_subsystems->suggestNext(time_current);
        if (!(time_next > time_current.getTime())) {
            //error
            ERR("error - the suggested next time step is invalid.");
            break;
        }
        TimeStep t_n1(time_current, time_next);
        INFO("\n#############################################################");
        INFO("Time step %d: t=%f s, dt=%f s ", t_n1.getTimeStepCount(), time_next, t_n1.getTimeStepSize());
        INFO("\n#############################################################");
        bool isAccepted = (_root_subsystems->solveTimeStep(t_n1)==0);
        if (isAccepted) {
            _root_subsystems->accept(time_next);
            time_current.accept(time_next);
            doSomethingAfterTimeStepAccepted(time_current);
        }
    }

    return time_current.getTimeStepCount();
}

}