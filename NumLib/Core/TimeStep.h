
#pragma once

namespace NumLib
{
typedef double TimeStep;
//class TimeStep;
class TimeStepping
{
public:
    void initialize();
    TimeStep next();
    void setNextTimeStep(TimeStep);
};

}