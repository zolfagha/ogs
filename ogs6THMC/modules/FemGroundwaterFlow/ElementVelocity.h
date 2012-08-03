/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ElementVelocity.h
 *
 * Created on 2012-07-13 by Norihiro Watanabe
 */

#pragma once

#include "BaseLib/CodingTools.h"
#include "DiscreteLib/Core/DiscreteSystem.h"
#include "FemLib/Function/FemFunction.h"
#include "NumLib/Function/Function.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "MaterialLib/PorousMedia.h"
#include "SolutionLib/FemProblem/AbstractTimeIndependentFemFunction.h"

#include "ProcessLib/TemplateTimeIndependentProcess.h"

class FunctionElementVelocity
    : public ProcessLib::TemplateTimeIndependentProcess<1,1>
{
    enum In { Head=0 };
    enum Out { Velocity=0 };
public:

    FunctionElementVelocity() 
    : _dis(0), _vel(0)
    {
    };

    virtual ~FunctionElementVelocity() {};


    virtual bool initialize(const BaseLib::Options &op);

    virtual void finalize() {};

    int solveTimeStep(const NumLib::TimeStep &/*time*/);

    virtual void accept(const NumLib::TimeStep &/*time*/);

private:
    DiscreteLib::DiscreteSystem* _dis;
    FemLib::FEMIntegrationPointFunctionVector* _vel;
    //NumLib::ITXFunction* _K;

    DISALLOW_COPY_AND_ASSIGN(FunctionElementVelocity);
};




