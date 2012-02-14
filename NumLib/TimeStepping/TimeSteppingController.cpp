
#include "TimeSteppingController.h"

#include <iostream>

namespace NumLib
{

void TimeSteppingController::setBeginning(double time) 
{
    _time_begin = time;
};

void TimeSteppingController::addTransientSystem(ITransientProblem &system) 
{
    _root_subsystems = &system;
};

void TimeSteppingController::solve(double time_end) 
{
    TimeStep time_current(_time_begin);
    std::cout << "\n\n***Start time steps\n";

    while (time_current.getTime()<time_end) {
        double time_next = _root_subsystems->suggestNext(time_current);
        TimeStep t_n1(time_current, time_next);
        bool isAccepted = (_root_subsystems->solve(t_n1)==0);
        if (isAccepted) {
            time_current.accept(time_next);
        }
    }
}

}