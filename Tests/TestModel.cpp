//
//#include <gtest/gtest.h>
//
//#include "ModelLib/GROUNDWATER_FLOW.h"
//
//using namespace ModelLib;
//
//TEST(Model, testGW)
//{
//    GROUNDWATER_FLOW gw;
//}

#include "NumLib/Core/TimeStep.h"

void test()
{
    // define problems and solution strategy
    // pass it to discretization systems
    // start time stepping
    double end_time = 1.;
    size_t time_step_cnt = 0;
    double current_time = 0;
    while ( current_time < end_time )
    {
        // seek next time step
        double next_time = .0;
        double dt = next_time - current_time; 
        // do something and check if this time step is ok
        bool accepted = false;
        if (accepted) {
            //
            current_time = next_time;
            ++time_step_cnt;
            // post process
        }
    }

}
