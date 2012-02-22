
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

void TimeSteppingController::solve(double time_end) 
{
    TimeStep time_current(_time_begin);
    std::cout << "\n\n***Start time steps\n";

    while (time_current.getTime()<time_end) {
        double time_next = _root_subsystems->suggestNext(time_current);
        if (!(time_next > time_current.getTime())) {
            //error
            std::cout << "error - the suggested next time step is invalid." << std::endl;
            break;
        }
        TimeStep t_n1(time_current, time_next);
        bool isAccepted = (_root_subsystems->solveTimeStep(t_n1)==0);
        if (isAccepted) {
            _root_subsystems->accept(time_next);
            time_current.accept(time_next);
        }
    }
}

}