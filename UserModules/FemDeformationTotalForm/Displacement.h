/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Displacement.h
 *
 * Created on 2012-07-13 by Norihiro Watanabe
 */

#pragma once

#include "BaseLib/CodingTools.h"
#include "MathLib/LinAlg/LinearEquations/LisInterface.h"
#include "MathLib/LinAlg/LinearEquations/DenseLinearEquations.h"
#include "FemLib/Function/FemIntegrationPointFunction.h"
#include "SolutionLib/Fem/FemEquation.h"
#include "SolutionLib/Fem/FemIVBVProblem.h"
#include "SolutionLib/Fem/SingleStepFEM.h"

#include "ProcessLib/TemplateTransientProcess.h"

#include "LinearElasticLinearLocalAssembler.h"
#include "LinearElasticResidualLocalAssembler.h"
#include "LinearElasticJacobianLocalAssembler.h"

#include "DiscreteLib/Serial/DiscreteSystem.h"
typedef DiscreteLib::DiscreteSystem MyDiscreteSystem;

typedef SolutionLib::TemplateFemEquation<
        MyDiscreteSystem,
        FemLinearElasticLinearLocalAssembler,
        FemLinearElasticResidualLocalAssembler,
        FemLinearElasticJacobianLocalAssembler
        >
        FemLinearElasticEquation;

typedef SolutionLib::FemIVBVProblem<
        MyDiscreteSystem,
        FemLinearElasticEquation
        > FemLinearElasticProblem;



class FunctionDisplacement
: public ProcessLib::TemplateTransientProcess<0,1>
{
    enum Out { Displacement=0 };
    //enum Out { Displacement=0, Strain=1, Stress=2 };
public:
    typedef FemLinearElasticProblem MyProblemType;
    typedef MyProblemType::EquationType MyEquationType;
    typedef SolutionLib::SingleStepFEM
            <
                FemLinearElasticProblem,
                //MathLib::DenseLinearEquations
                MathLib::CRSLisSolver
            > MySolutionType;
    typedef FemLib::FemNodalFunctionVector<MyDiscreteSystem>::type MyNodalFunctionVector;

    ///
    FunctionDisplacement()
    : ProcessLib::TemplateTransientProcess<0,1>("DEFORMATION"),
      _problem(0), _solution(0), _feObjects(0), _displacement(0)//, _strain(0), _stress(0)
    {
    };

    ///
    virtual ~FunctionDisplacement()
    {
        BaseLib::releaseObject(_problem, _solution, _feObjects, _displacement /*, _strain, _stress*/);
    }

    /// initialize this process
    virtual bool initialize(const BaseLib::Options &option);

    /// finalize but nothing to do here
    virtual void finalize() {};


protected:
    virtual void updateOutputParameter(const NumLib::TimeStep &time);

    virtual MySolutionType* getSolution() {return _solution;};

    virtual void output(const NumLib::TimeStep &time);

private:
    //void calculateStressStrain();

    DISALLOW_COPY_AND_ASSIGN(FunctionDisplacement);

private:
    FemLinearElasticProblem* _problem;
    MySolutionType* _solution;
    FemLib::LagrangianFeObjectContainer* _feObjects;
    MyNodalFunctionVector* _displacement;
//    FemLib::FEMIntegrationPointFunctionVector* _strain;
//    FemLib::FEMIntegrationPointFunctionVector* _stress;

};



