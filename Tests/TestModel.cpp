
#include <gtest/gtest.h>

//#include "ModelLib/GROUNDWATER_FLOW.h"
//
//using namespace ModelLib;
//
//TEST(Model, testGW)
//{
//    GROUNDWATER_FLOW gw;
//}

#include "NumLib/TimeStepping/TimeSteppingController.h"

using namespace NumLib;


TEST(Model, test1)
{
    // define problems and solution strategy
    ITransientSystem* problem;

    // pass it to discretization systems
    TimeSteppingController timeStepping;
    timeStepping.addTransientSystem(*problem);

    // start time stepping
    timeStepping.setBeginning(.0);
    timeStepping.solve(100.);
}
