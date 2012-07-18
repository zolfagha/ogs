
#pragma once

#include <string>
#include <vector>

#include "IOutputTiming.h"

class OutputTimingBuilder
{
public:
    IOutputTiming* create(const std::string &name, size_t n, std::vector<double> *vec_time = 0);
};
