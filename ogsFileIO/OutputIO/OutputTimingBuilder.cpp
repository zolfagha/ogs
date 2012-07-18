
#include "OutputTimingBuilder.h"

#include "OutputTimingStepPeriodic.h"

IOutputTiming* OutputTimingBuilder::create(const std::string &name, size_t n, std::vector<double> *vec_time)
{
    if (name.compare("STEPS")==0) {
        return new OutputTimingStepPeriodic(n);
    }
    return NULL;
}
