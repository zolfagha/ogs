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
#include "DiscreteLib/Serial/DiscreteSystem.h"
#include "FemLib/Function/FemFunction.h"
#include "NumLib/TimeStepping/TimeStep.h"
#include "ProcessLib/TemplateTimeIndependentProcess.h"

typedef DiscreteLib::DiscreteSystem MyDiscreteSystem;

/**
 * \brief Stress, strain evaluator based on total displacements
 */
class FunctionElementStressStrain
    : public ProcessLib::TemplateTimeIndependentProcess<1,2>
{
    enum In { Displacement=0 };
    enum Out { Strain=0, Stress=1};
    typedef FemLib::FemNodalFunctionVector<MyDiscreteSystem>::type MyNodalFunctionVector;
    typedef FemLib::FEMIntegrationPointFunctionVector<MyDiscreteSystem>::type MyIntegrationPointFunctionVector;
    typedef NumLib::TXVectorFunctionAsColumnData<MyIntegrationPointFunctionVector> IntegrationPointScalarWrapper;
    typedef NumLib::TXVectorFunctionAsColumnData<MyNodalFunctionVector> NodalPointScalarWrapper;
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
    MyDiscreteSystem* _dis;
    MyIntegrationPointFunctionVector* _strain;
    MyIntegrationPointFunctionVector* _stress;
    std::vector<IntegrationPointScalarWrapper*> _vec_strain_components;
    std::vector<IntegrationPointScalarWrapper*> _vec_stress_components;
};
