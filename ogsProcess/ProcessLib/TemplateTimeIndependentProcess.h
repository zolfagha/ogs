
#pragma once

#include "BaseLib/Options.h"
#include "NumLib/TransientCoupling/TransientMonolithicSystem.h"
#include "Process.h"

namespace ProcessLib
{

template <
	size_t N_IN_PARAMETER,
	size_t N_OUT_PARAMETER
	>
class TemplateTimeIndependentProcess
: public Process
{
public:
	TemplateTimeIndependentProcess()
    {
        AbstractTransientMonolithicSystem::resizeInputParameter(N_IN_PARAMETER);
        AbstractTransientMonolithicSystem::resizeOutputParameter(N_OUT_PARAMETER);
    }

    double suggestNext(const NumLib::TimeStep &/*time_current*/) { return .0; }

    bool isAwake(const NumLib::TimeStep &/*time*/) { return true;  }

    virtual void accept(const NumLib::TimeStep &/*time*/) {};

};
}

