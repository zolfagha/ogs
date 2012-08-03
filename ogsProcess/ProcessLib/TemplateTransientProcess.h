/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file TemplateTransientProcess.h
 *
 * Created on 2012-08-03 by Norihiro Watanabe
 */

#pragma once

#include "logog/include/logog.hpp"

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
    explicit TemplateTransientProcess(const std::string &pcs_name) 
        : _pcs_name(pcs_name)
    {
        Process::resizeInputParameter(N_IN_PARAMETER);
        Process::resizeOutputParameter(N_OUT_PARAMETER);
    }

    ///
    virtual ~TemplateTransientProcess() {};

    ///
    virtual int solveTimeStep(const NumLib::TimeStep &time)
    {
        INFO("Solving %s...", _pcs_name.c_str());
        initializeTimeStep(time);
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
    virtual void initializeTimeStep(const NumLib::TimeStep &/*time*/) {};
    ///
    virtual void updateOutputParameter(const NumLib::TimeStep &time) = 0;
    ///
    virtual void output(const NumLib::TimeStep &time) = 0;

private:
    std::string _pcs_name;
};

} //end ProcessLib

