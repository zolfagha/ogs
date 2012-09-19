/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Concentration.h
 *
 * Created on 2012-09-06 by Haibing Shao
 */

#ifndef CONCENTRATIONS_H
#define CONCENTRATIONS_H

#include "BaseLib/CodingTools.h"
#include "NumLib/TransientAssembler/ElementWiseTimeEulerEQSLocalAssembler.h"
#include "NumLib/TransientAssembler/ElementWiseTimeEulerResidualLocalAssembler.h"
#include "NumLib/Function/DiscreteDataConvergenceCheck.h"
#include "SolutionLib/Fem/FemIVBVProblem.h"
#include "SolutionLib/Fem/TemplateFemEquation.h"
#include "SolutionLib/Fem/SingleStepFEM.h"
#include "ProcessLib/AbstractTransientProcess.h"


template <class T_DISCRETE_SYSTEM, class T_LINEAR_SOLVER>
class FunctionConcentrations
	: public ProcessLib::Process
{
public:
    enum In { Velocity=0 };
    enum Out { Concentrations = 0 };

    typedef T_DISCRETE_SYSTEM MyDiscreteSystem;
    typedef T_LINEAR_SOLVER MyLinearSolver;

    // memory for discretized concentration vector
    typedef typename FemLib::FemNodalFunctionVector<MyDiscreteSystem>::type MyNodalFunctionVector;

    // local assembler
	// for the linear systems, use the same settings as Mass_Transport
    typedef MassTransportTimeODELocalAssembler<NumLib::ElementWiseTimeEulerEQSLocalAssembler> MyLinearAssemblerType;          // TO BE CHANGED
    typedef MassTransportTimeODELocalAssembler<NumLib::ElementWiseTimeEulerResidualLocalAssembler> MyResidualAssemblerType;   // TO BE CHANGED
    typedef MassTransportJacobianLocalAssembler MyJacobianAssemblerType;                                                      // TO BE CHANGED
	// for the nonlinear part, use different settings
	
	// Equation definition
    typedef SolutionLib::TemplateFemEquation<
            MyDiscreteSystem,
            MyLinearSolver,
            MyLinearAssemblerType,
            MyResidualAssemblerType,
            MyJacobianAssemblerType
            > MyLinearEquationType;
    // IVBV problem definition
    typedef SolutionLib::FemIVBVProblem<
            MyDiscreteSystem,
            MyLinearEquationType
            > MyLinearTransportProblemType;
    // Solution algorithm definition
    typedef SolutionLib::SingleStepFEM<
            MyLinearTransportProblemType,
            MyLinearSolver
            > MyLinearSolutionType;

	// nonlinear coupled part
	//// Equation definition
	//   typedef SolutionLib::TemplateFemEquation<
	//           MyDiscreteSystem,
	//           MyLinearSolver,
	//           MyLinearAssemblerType,
	//           MyResidualAssemblerType,
	//           MyJacobianAssemblerType
	//           > MyLinearEquationType;
	//   // IVBV problem definition
	//   typedef SolutionLib::FemIVBVProblem<
	//           MyDiscreteSystem,
	//           MyLinearEquationType
	//           > MyLinearTransportProblemType;
	//   // Solution algorithm definition
	//   typedef SolutionLib::SingleStepFEM<
	//           MyLinearTransportProblemType,
	//           MyLinearSolver
	//           > MyLinearSolutionType;

    FunctionConcentrations() 
        : Process("KIN_REACT_GIA", 1, 1),
          _feObjects(0)
    {
        // set default parameter name
		ProcessLib::Process::setInputParameterName(Velocity, "Velocity");
        ProcessLib::Process::setOutputParameterName(Concentrations, "Concentrations");
    };

    virtual ~FunctionConcentrations()
    {
        BaseLib::releaseObject(_feObjects);
        BaseLib::releaseObject(_linear_problem, _linear_solution); 
        BaseLib::releaseObject(_ReductionKin);
        BaseLib::releaseObject(_concentrations, _eta, _xi);
    };

    /// initialize this process
    virtual bool initialize(const BaseLib::Options &option);

    /// finalize but nothing to do here
    virtual void finalize() {};

    virtual NumLib::IConvergenceCheck* getConvergenceChecker()
    {
        return &_checker;
    }

	virtual int solveTimeStep(const NumLib::TimeStep &time)
    {
		// TODO
        //INFO("Solving %s...", getProcessName().c_str());
        //initializeTimeStep(time);
        //getSolution()->solveTimeStep(time);
        //updateOutputParameter(time);
        return 0;
    }

    /// 
    virtual double suggestNext(const NumLib::TimeStep &time_current) 
    {
		// TODO: only one solution is enough. 
        return getSolution()->suggestNext(time_current); 
    }

    ///
    virtual bool isAwake(const NumLib::TimeStep &time) 
    { 
		// TODO: only one solution is enough
        return getSolution()->isAwake(time);  
    }

    ///
    virtual void accept(const NumLib::TimeStep &time)
    {
        output(time);
		// TODO: call accept for all solutions. 
        // getSolution()->accept(time);
    };

protected:
    virtual void initializeTimeStep(const NumLib::TimeStep &time);

    virtual void updateOutputParameter(const NumLib::TimeStep &time);

    virtual MyLinearSolutionType* getSolution() {return _linear_solution;};

    virtual void output(const NumLib::TimeStep &time);

private:
    DISALLOW_COPY_AND_ASSIGN(FunctionConcentrations);

private:
    // linear problem
    MyLinearTransportProblemType* _linear_problem;
    // linear solver
    MyLinearSolutionType* _linear_solution;

    FemLib::LagrangianFeObjectContainer* _feObjects; 
    
	NumLib::DiscreteDataConvergenceCheck _checker;

    // pointer to the reduction scheme. 
	ogsChem::chemReductionKin* _ReductionKin; 

    // concentrations vector, including all components in the MCP data structure
    MyNodalFunctionVector* _concentrations; 
    // eta vector, including eta_mobile and eta_immobile
    MyNodalFunctionVector* _eta;
    // xi vector, including xi_mobile and xi_immobile parts
    MyNodalFunctionVector* _xi; 
}; 

#include "Concentrations.hpp"

#endif  // end of ifndef