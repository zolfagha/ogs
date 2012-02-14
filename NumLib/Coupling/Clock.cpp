
#include "Clock.h"

namespace NumLib
{

void TimeSteppingController::setBeginning(TimeStep time) 
{
    _time_begin = time;
};

void TimeSteppingController::addTransientSystem(ITransientCoupledProblem &system) 
{
    _root_subsystems = &system;
};

void TimeSteppingController::solve(TimeStep time_end) 
{
    TimeStep time_current = _time_begin;
    while (time_current<time_end) {
        TimeStep try_next = _root_subsystems->suggestNext(time_current);
        while (_root_subsystems->solve(try_next)!=0) {
            try_next = _root_subsystems->suggestNext(time_current);
        }
        time_current = try_next;
    }
}

}