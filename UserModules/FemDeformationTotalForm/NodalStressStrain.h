/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file NodalStressStrain.h
 *
 * Created on 2012-08-14 by Norihiro Watanabe
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
class FunctionNodalStressStrain
    : public ProcessLib::TemplateTimeIndependentProcess<2,2>
{
    enum In { GpStrain=0, GpStress=1};
    enum Out { NodStrain=0, NodStress=1};
    typedef FemLib::FemNodalFunctionVector<MyDiscreteSystem>::type MyNodalFunctionVector;
    typedef FemLib::FEMIntegrationPointFunctionVector<MyDiscreteSystem>::type MyIntegrationPointFunctionVector;
    typedef NumLib::TXVectorFunctionAsColumnData<MyIntegrationPointFunctionVector> IntegrationPointScalarWrapper;
    typedef NumLib::TXVectorFunctionAsColumnData<MyNodalFunctionVector> NodalPointScalarWrapper;
public:

    FunctionNodalStressStrain() : _dis(0), _nodal_strain(0), _nodal_stress(0)
    {
    };

    virtual ~FunctionNodalStressStrain()
    {
        BaseLib::releaseObject(_nodal_strain, _nodal_stress);
        BaseLib::releaseObjectsInStdVector(_vec_nodal_strain_components);
        BaseLib::releaseObjectsInStdVector(_vec_nodal_stress_components);
    };

    virtual bool initialize(const BaseLib::Options &op);

    virtual void finalize() {};

    int solveTimeStep(const NumLib::TimeStep &/*time*/);

    void accept(const NumLib::TimeStep &/*time*/);

private:
    size_t getNumberOfStrainComponents() const;

    DISALLOW_COPY_AND_ASSIGN(FunctionNodalStressStrain);

private:
    MyDiscreteSystem* _dis;
    MyNodalFunctionVector* _nodal_strain;
    MyNodalFunctionVector* _nodal_stress;
    std::vector<NodalPointScalarWrapper*> _vec_nodal_strain_components;
    std::vector<NodalPointScalarWrapper*> _vec_nodal_stress_components;
};
