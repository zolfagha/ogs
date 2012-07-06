
#pragma once

#include "NumLib/TransientCoupling/TransientMonolithicSystem.h"
#include "SolutionLib/Solution/SingleStepFEM.h"

namespace SolutionLib
{

template <
	size_t N_IN_PARAMETER,
	size_t N_OUT_PARAMETER
	>
class AbstractTransientFemFunction
: public NumLib::TemplateTransientMonolithicSystem
{
public:
    AbstractTransientFemFunction()
    {
        TemplateTransientMonolithicSystem::resizeInputParameter(N_IN_PARAMETER);
        TemplateTransientMonolithicSystem::resizeOutputParameter(N_OUT_PARAMETER);
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
    virtual AbstractTimeSteppingAlgorithm* getSolution() = 0;
    virtual void updateOutput() = 0;
};

}
