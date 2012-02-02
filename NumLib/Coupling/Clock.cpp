
#include "Clock.h"

namespace NumLib
{
void Clock::setBeginning(TimeStep time) 
{
    _time_begin = time;
};

void Clock::addTransientSystem(ITransientSystem *system) 
{
    rootAsyncPartSolution.addChildren(system);
};

void Clock::moveForwardUntill(TimeStep time_end) 
{
    if (rootAsyncPartSolution.getNumberOfChildren()==0)
        return;

    TimeStep time_current = _time_begin;
    while (time_current<time_end) {
        TimeStep try_next = rootAsyncPartSolution.suggestNext(time_current);
        while (!rootAsyncPartSolution.solveNextStep(try_next)) {
            try_next = rootAsyncPartSolution.suggestNext(time_current);
        }
        time_current = try_next;
    }
}

}