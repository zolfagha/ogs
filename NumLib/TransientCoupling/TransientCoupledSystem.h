
#pragma once

#include "NumLib/TimeStepping/TimeStep.h"
#include "NumLib/TimeStepping/ITransientSystem.h"
#include "NumLib/Coupling/ICoupledProblem.h"

namespace NumLib
{

/**
 * \brief Interface class of transient coupling problems
 */
class ITransientCoupledSystem : public ICoupledSystem, public ITransientSystem
{
public:
    int solve(const TimeStep &time) 
    {
        setCurrentTime(time);
        return solve();
    }
private:
    int solve()
    {
        return solveTimeStep(getCurrentTime());
    }
};

}
