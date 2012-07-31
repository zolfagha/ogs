
#pragma once

#include <string>
#include <vector>
#include <map>

#include "BaseLib/CodingTools.h"
#include "NumLib/Function/TXFunction.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "CommonIO/PVDWriter.h"

#include "IOutput.h"

class PVDOutput : public IOutput
{
public:
    PVDOutput() : _pvd(0) {};
    virtual ~PVDOutput() 
    {
        BaseLib::releaseObject(_pvd);
    }
    virtual void write( const NumLib::TimeStep &current_time,
                        BaseLib::OrderedMap<std::string, OutputVariableInfo> &data);

private:
    PVDWriter* _pvd;
};
