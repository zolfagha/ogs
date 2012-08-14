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

#include <vector>

#include "BaseLib/CodingTools.h"
#include "BaseLib/Options.h"
#include "DiscreteLib/Core/DiscreteSystem.h"
#include "FemLib/Function/FemFunction.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "ProcessLib/TemplateTimeIndependentProcess.h"

/**
 * \brief Stress, strain evaluator based on total displacements
 */
class FunctionElementStressStrain
    : public ProcessLib::TemplateTimeIndependentProcess<1,2>
{
    enum In { Displacement=0 };
    enum Out { Strain=0, Stress=1};
    typedef NumLib::TXVectorFunctionAsColumnData<FemLib::FEMIntegrationPointFunctionVector> IntegrationPointScalarWrapper;
    typedef NumLib::TXVectorFunctionAsColumnData<FemLib::FemNodalFunctionVector> NodalPointScalarWrapper;
public:

    FunctionElementStressStrain() : _dis(0), _strain(0), _stress(0)
    {
    };

    virtual ~FunctionElementStressStrain()
    {
        BaseLib::releaseObject(_strain, _stress);
        BaseLib::releaseObjectsInStdVector(_vec_strain_components);
        BaseLib::releaseObjectsInStdVector(_vec_stress_components);
    };

    virtual bool initialize(const BaseLib::Options &op);

    virtual void finalize() {};

    int solveTimeStep(const NumLib::TimeStep &/*time*/);

    void accept(const NumLib::TimeStep &/*time*/);

private:
    size_t getNumberOfStrainComponents() const;

    DISALLOW_COPY_AND_ASSIGN(FunctionElementStressStrain);

private:
    DiscreteLib::DiscreteSystem* _dis;
    FemLib::FEMIntegrationPointFunctionVector* _strain;
    FemLib::FEMIntegrationPointFunctionVector* _stress;
    std::vector<IntegrationPointScalarWrapper*> _vec_strain_components;
    std::vector<IntegrationPointScalarWrapper*> _vec_stress_components;
};
