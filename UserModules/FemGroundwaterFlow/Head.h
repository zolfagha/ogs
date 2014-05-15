/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Head.h
 *
 * Created on 2012-07-13 by Norihiro Watanabe
 */

#pragma once

#include "BaseLib/CodingTools.h"
#include "NumLib/Function/DiscreteDataConvergenceCheck.h"
#include "NumLib/TransientAssembler/ElementWiseTimeEulerEQSLocalAssembler.h"
#include "NumLib/TransientAssembler/ElementWiseTimeEulerResidualLocalAssembler.h"
#include "SolutionLib/Fem/FemIVBVProblem.h"
#include "SolutionLib/Fem/TemplateFemEquation.h"
#include "SolutionLib/Fem/SingleStepFEM.h"
#include "ProcessLib/AbstractTransientProcess.h"

#include "HeadBasedGWTimeODELocalAssembler.h"
#include "HeadBasedGWJacobianLocalAssembler.h"


/**
 * \brief Hydraulic head calculator based on Groundwater flow equation using FEM
 *
 * Output result
 * - Hydraulic head (nodal)
 *
 * \tparam T_DISCRETE_SYSTEM    Discrete system type
 * \tparam T_LINEAR_SOLVER      Linear solver type
 */
template <class T_DISCRETE_SYSTEM, class T_LINEAR_SOLVER>
class FunctionHead
: public ProcessLib::AbstractTransientProcess
{
public:
    // define input and output parameter id here
    enum Out { Head=0 };

    //
    typedef T_DISCRETE_SYSTEM MyDiscreteSystem;
    typedef T_LINEAR_SOLVER MyLinearSolver;
    // local assembler
    typedef HeadBasedGWTimeODELocalAssembler<
                NumLib::ElementWiseTimeEulerEQSLocalAssembler
                > MyLinearAssemblerType;
    typedef HeadBasedGWTimeODELocalAssembler<
                NumLib::ElementWiseTimeEulerResidualLocalAssembler
                > MyResidualAssemblerType;
    typedef HeadBasedGWJacobianLocalAssembler MyJacobianAssemblerType;
    // Equation definition
    typedef SolutionLib::TemplateFemEquation<
            MyDiscreteSystem,
            MyLinearSolver,
            MyLinearAssemblerType,
            MyResidualAssemblerType,
            MyJacobianAssemblerType
            >
            MyEquationType;
    // IVBV problem definition
    typedef SolutionLib::FemIVBVProblem<
            MyDiscreteSystem,
            MyEquationType
            > MyProblemType;
    // Solution algorithm definition
    typedef SolutionLib::SingleStepFEM
            <
            MyProblemType,
            MyLinearSolver
            > MySolutionType;

    typedef typename FemLib::FemNodalFunctionScalar<MyDiscreteSystem>::type MyNodalFunctionScalar;
    typedef typename MyProblemType::MyVariable MyVariable;
    
    ///
    FunctionHead() 
    : AbstractTransientProcess("GROUNDWATER_FLOW", 0, 1), _problem(0), _solution(0), _feObjects(0)
    {
        // set default parameter name
        ProcessLib::AbstractTransientProcess::setOutputParameterName(Head, "Head");
        _checker.addValueType(NumLib::DiscreteDataType::NodalScalar);
    };

    ///
    virtual ~FunctionHead()
    {
        BaseLib::releaseObject(_problem, _solution, _feObjects);
    }

    /// initialize this process
    virtual bool initialize(const BaseLib::Options &option);

    /// finalize but nothing to do here
    virtual void finalize() {};

    /// return a convergence check class for this function
    virtual NumLib::IConvergenceCheck* getConvergenceChecker()
    {
        return &_checker;
    }

protected:
    virtual void updateOutputParameter(const NumLib::TimeStep &time);

    virtual MySolutionType* getSolution() {return _solution;};

    virtual void output(const NumLib::TimeStep &time);


private:
    DISALLOW_COPY_AND_ASSIGN(FunctionHead);

private:
    MyProblemType* _problem;
    MySolutionType* _solution;
    FemLib::LagrangeFeObjectContainer* _feObjects;
    NumLib::DiscreteDataConvergenceCheck _checker;
};

#include "Head.hpp"


