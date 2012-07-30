
#pragma once

#include "BaseLib/Options.h"
#include "NumLib/TransientCoupling/TransientMonolithicSystem.h"
#include "SolutionLib/Solution/AbstractTimeSteppingAlgorithm.h"
#include "Process.h"

namespace ProcessLib
{

/**
 * \brief Implementation of Process (ITransientSystem) class for monolithic system
 *
 * \tparam N_IN_PARAMETER the number of input parameters 
 * \tparam N_OUT_PARAMETER the number of output parameters 
 */
template <
    size_t N_IN_PARAMETER,
    size_t N_OUT_PARAMETER
    >
class TemplateTransientProcess : public Process
{
public:
    ///
    TemplateTransientProcess()
    {
        AbstractTransientMonolithicSystem::resizeInputParameter(N_IN_PARAMETER);
        AbstractTransientMonolithicSystem::resizeOutputParameter(N_OUT_PARAMETER);
    }

    ///
    virtual ~TemplateTransientProcess() {};

    ///
    virtual int solveTimeStep(const NumLib::TimeStep &time)
    {
        getSolution()->solveTimeStep(time);
        updateOutputParameter(time);
        return 0;
    }

    /// 
    virtual double suggestNext(const NumLib::TimeStep &time_current) 
    {
        return getSolution()->suggestNext(time_current); 
    }

    ///
    virtual bool isAwake(const NumLib::TimeStep &time) 
    { 
        return getSolution()->isAwake(time);  
    }

    ///
    virtual void accept(const NumLib::TimeStep &time)
    {
        output(time);
        getSolution()->accept(time);
    };

protected:
    ///
    virtual SolutionLib::AbstractTimeSteppingAlgorithm* getSolution() = 0;
    ///
    virtual void updateOutputParameter(const NumLib::TimeStep &time) = 0;
    ///
    virtual void output(const NumLib::TimeStep &time) = 0;
};

} //end ProcessLib

