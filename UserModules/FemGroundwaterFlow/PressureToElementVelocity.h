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
#include "NumLib/Function/DiscreteDataConvergenceCheck.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "ProcessLib/AbstractTimeIndependentProcess.h"

template <class T_DISCRETE_SYSTEM>
class FunctionPressureToElementVelocity
    : public ProcessLib::AbstractTimeIndependentProcess
{
public:
    enum In { Pressure=0 };
    enum Out { Velocity=0 };

    typedef T_DISCRETE_SYSTEM MyDiscreteSystem;
    typedef typename FemLib::FemNodalFunctionScalar<MyDiscreteSystem>::type MyNodalFunctionScalar;
    typedef typename FemLib::FEMIntegrationPointFunctionVector<MyDiscreteSystem>::type MyIntegrationPointFunctionVector;

    FunctionPressureToElementVelocity()
    : ProcessLib::AbstractTimeIndependentProcess("PRESSURE_TO_ELEMENT_VELOCITY", 1, 1), _dis(NULL), _vel(NULL), _feObjects(NULL)
    {
        // set default parameter name
        ProcessLib::AbstractTimeIndependentProcess::setInputParameterName(Pressure, "Pressure");
        ProcessLib::AbstractTimeIndependentProcess::setOutputParameterName(Velocity, "Velocity");
    };

    virtual ~FunctionPressureToElementVelocity()
    {
        BaseLib::releaseObject(_feObjects, _vel);
    };


    virtual bool initialize(const BaseLib::Options &op);

    virtual void finalize() {};

    int solveTimeStep(const NumLib::TimeStep &/*time*/);

    virtual void accept(const NumLib::TimeStep &/*time*/);

    ///
    virtual NumLib::IConvergenceCheck* getConvergenceChecker() { return &_checker; };

private:
    DiscreteLib::IDiscreteSystem* _dis;
    MyIntegrationPointFunctionVector* _vel;
    FemLib::LagrangianFeObjectContainer* _feObjects;
    NumLib::DiscreteDataConvergenceCheck _checker;

    DISALLOW_COPY_AND_ASSIGN(FunctionPressureToElementVelocity);
};

#include "PressureToElementVelocity.hpp"



