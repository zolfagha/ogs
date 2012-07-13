
#pragma once

#include "BaseLib/Options.h"
#include "NumLib/TransientCoupling/TransientMonolithicSystem.h"
#include "SolutionLib/Solution/AbstractTimeSteppingAlgorithm.h"
#include "Process.h"

namespace ProcessLib
{

template <
	size_t N_IN_PARAMETER,
	size_t N_OUT_PARAMETER
	>
class TemplateTransientProcess : public Process
{
public:
	TemplateTransientProcess()
    {
		AbstractTransientMonolithicSystem::resizeInputParameter(N_IN_PARAMETER);
		AbstractTransientMonolithicSystem::resizeOutputParameter(N_OUT_PARAMETER);
    }

    int solveTimeStep(const NumLib::TimeStep &time)
    {
    	getSolution()->solveTimeStep(time);
        updateOutput();
        return 0;
    }

    double suggestNext(const NumLib::TimeStep &time_current) { return getSolution()->suggestNext(time_current); }

    bool isAwake(const NumLib::TimeStep &time) { return getSolution()->isAwake(time);  }

    void accept(const NumLib::TimeStep &time)
    {
    	getSolution()->accept(time);
    };

protected:
    virtual SolutionLib::AbstractTimeSteppingAlgorithm* getSolution() = 0;
    virtual void updateOutput() = 0;
};

}

