
#include <gtest/gtest.h>

#include <vector>

#include "NumLib/Coupling/Clock.h"
#include "NumLib/Coupling/TransientSystems.h"

using namespace NumLib;

TEST(Num, test1)
{
    std::vector<ITransientSystem*> vec_systems;

    Clock clock;
    for (size_t i=0; i<vec_systems.size(); i++)
        clock.addTransientSystem(vec_systems[i]);

    TimeStep t0, t_end;
    clock.setBeginning(t0);

    // start clock
    clock.moveForwardUntill(t_end); 
}

