/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file StressStrain.h
 *
 * Created on 2012-07-13 by Norihiro Watanabe
 */

#pragma once

#include "BaseLib/CodingTools.h"
#include "BaseLib/Options.h"
#include "DiscreteLib/Core/DiscreteSystem.h"
#include "FemLib/Function/FemFunction.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "ProcessLib/TemplateTimeIndependentProcess.h"

/**
 * \brief Stress, strain evaluator based on total displacements
 */
class FunctionStressStrain
    : public ProcessLib::TemplateTimeIndependentProcess<1,2>
{
    enum In { Displacement=0 };
    enum Out { Stress=0, Strain=1};

public:

    FunctionStressStrain() : _dis(0), _strain(0), _stress(0)
    {
    };

    virtual ~FunctionStressStrain() {};

    virtual bool initialize(const BaseLib::Options &op);

    virtual void finalize() {};

    int solveTimeStep(const NumLib::TimeStep &/*time*/);

    void accept(const NumLib::TimeStep &/*time*/);

private:
    size_t getNumberOfStrainComponents() const;

    DISALLOW_COPY_AND_ASSIGN(FunctionStressStrain);

private:
    DiscreteLib::DiscreteSystem* _dis;
    FemLib::FEMIntegrationPointFunctionVector* _strain;
    FemLib::FEMIntegrationPointFunctionVector* _stress;
};
