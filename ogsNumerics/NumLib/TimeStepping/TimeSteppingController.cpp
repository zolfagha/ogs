
#include "TimeSteppingController.h"

#include <iostream>

namespace NumLib
{

void TimeSteppingController::setBeginning(double time) 
{
    _time_begin = time;
};

void TimeSteppingController::addTransientSystem(ITransientSystem &system) 
{
    _root_subsystems = &system;
};

size_t TimeSteppingController::solve(double time_end) 
{
    TimeStep time_current(_time_begin);
    std::cout << "\n\n***Start time steps with t0=" << _time_begin << "s\n";

    while (time_current.getTime()<time_end) {
        double time_next = _root_subsystems->suggestNext(time_current);
        if (!(time_next > time_current.getTime())) {
            //error
            std::cout << "error - the suggested next time step is invalid." << std::endl;
            break;
        }
        TimeStep t_n1(time_current, time_next);
        std::cout << "\n#############################################################" << std::endl;
        std::cout << "Time step " << t_n1.getTimeStepID() << ": t=" << time_next <<  "s, dt=" << t_n1.getTimeStepSize() << "s " << std::endl  << std::endl;
        bool isAccepted = (_root_subsystems->solveTimeStep(t_n1)==0);
        if (isAccepted) {
            _root_subsystems->accept(time_next);
            time_current.accept(time_next);
        }
        std::cout << "#############################################################\n" << std::endl;
    }

    return time_current.getTimeStepID();
}

}