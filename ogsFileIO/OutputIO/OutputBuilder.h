
#pragma once

#include <string>

#include "IOutput.h"

class OutputBuilder
{
public:
    IOutput* create(const std::string &name);
};
