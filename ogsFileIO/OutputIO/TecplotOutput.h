
#pragma once

#include <string>
#include <vector>
#include <map>

#include "BaseLib/CodingTools.h"
#include "NumLib/Function/TXFunction.h"
#include "NumLib/TimeStepping/TimeStep.h"

#include "IOutput.h"

class TecplotOutput : public IOutput
{
public:
    virtual void write(const NumLib::TimeStep &current_time,BaseLib::OrderedMap<std::string, OutputVariableInfo> &data)
    {

    }
};
