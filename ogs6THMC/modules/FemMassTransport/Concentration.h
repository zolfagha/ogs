/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Concentration.h
 *
 * Created on 2012-07-13 by Norihiro Watanabe
 */

#pragma once

#include "BaseLib/CodingTools.h"
#include "MathLib/LinAlg/LinearEquations/LisInterface.h"
#include "NumLib/TransientAssembler/ElementWiseTimeEulerEQSLocalAssembler.h"
#include "NumLib/TransientAssembler/ElementWiseTimeEulerResidualLocalAssembler.h"
#include "SolutionLib/FemProblem/FemIVBVProblem.h"
#include "SolutionLib/FemProblem/AbstractTransientFemFunction.h"
#include "ProcessLib/TemplateTransientProcess.h"
#include "MaterialLib/Compound.h"

#include "MassTransportTimeODELocalAssembler.h"
#include "MassTransportJacobianLocalAssembler.h"

//--------------------------------------------------------------------------------------------------
// Equation definition
//--------------------------------------------------------------------------------------------------
typedef SolutionLib::TemplateFemEquation<
        MassTransportTimeODELocalAssembler<NumLib::ElementWiseTimeEulerEQSLocalAssembler>,
        MassTransportTimeODELocalAssembler<NumLib::ElementWiseTimeEulerResidualLocalAssembler>,
        MassTransportJacobianLocalAssembler
        >
        FemMTEquation;

//--------------------------------------------------------------------------------------------------
// IVBV problem definition
//--------------------------------------------------------------------------------------------------
typedef SolutionLib::FemIVBVProblem< FemMTEquation > FemMTProblem;


//--------------------------------------------------------------------------------------------------
// Function definition
//--------------------------------------------------------------------------------------------------
class FunctionConcentration
: public ProcessLib::TemplateTransientProcess<1,1>
{
    enum In { Velocity=0 };
    enum Out { Concentration = 0 };

public:
    typedef FemMTProblem MyProblemType;
    typedef MyProblemType::EquationType MyEquationType;
    typedef SolutionLib::SingleStepFEM
            <
                FemMTProblem,
                MathLib::CRSLisSolver
            > MySolutionType;

    FunctionConcentration() 
        : TemplateTransientProcess<1,1>("MASS_TRANSPORT"),
          _problem(0), _solution(0), _feObjects(0), _compound(0)
    {
    };

    virtual ~FunctionConcentration()
    {
        BaseLib::releaseObject(_problem, _solution, _feObjects);
    };

    /// initialize this process
    virtual bool initialize(const BaseLib::Options &option);

    /// finalize but nothing to do here
    virtual void finalize() {};

protected:
    virtual void initializeTimeStep(const NumLib::TimeStep &time);

    virtual void updateOutputParameter(const NumLib::TimeStep &time);

    virtual MySolutionType* getSolution() {return _solution;};

    virtual void output(const NumLib::TimeStep &time);

private:
    DISALLOW_COPY_AND_ASSIGN(FunctionConcentration);

private:
    MyProblemType* _problem;
    MySolutionType* _solution;
    FemLib::LagrangianFeObjectContainer* _feObjects;
    MaterialLib::Compound* _compound;
};

