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
#include "DiscreteLib/Core/IDiscreteSystem.h"
#include "FemLib/Function/FemFunction.h"
#include "NumLib/Function/Function.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "ProcessLib/TemplateTimeIndependentProcess.h"

#include "DiscreteLib/Serial/DiscreteSystem.h"
typedef DiscreteLib::DiscreteSystem MyDiscreteSystem;

class FunctionElementVelocity
    : public ProcessLib::TemplateTimeIndependentProcess<1,1>
{
    enum In { Head=0 };
    enum Out { Velocity=0 };
    typedef FemLib::FemNodalFunctionScalar<MyDiscreteSystem>::type MyNodalFunctionScalar;
    typedef FemLib::FEMIntegrationPointFunctionVector<MyDiscreteSystem>::type MyIntegrationPointFunctionVector;
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
    DiscreteLib::IDiscreteSystem* _dis;
    MyIntegrationPointFunctionVector* _vel;

    DISALLOW_COPY_AND_ASSIGN(FunctionElementVelocity);
};




