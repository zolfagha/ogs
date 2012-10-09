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
#include "FemKinReduction.h"
#include "SingleStepKinReduction.h"
#include "SolutionLib/Fem/TemplateFemEquation.h"
#include "SolutionLib/Fem/SingleStepFEM.h"
#include "ProcessLib/AbstractTransientProcess.h"
#include "LinearTransportTimeODELocalAssember.h"
#include "LinearTransportJacobianLocalAssembler.h"
#include "MathLib/DataType.h"
#include "NonLinearReactiveTransportTimeODELocalAssembler.h"
#include "NonLinearReactiveTransportJacabianLocalAssembler.h"
#include "NestedOdeNRIterationStepInitializer.h"
#include "NumLib/Nonlinear/DiscreteNRSolverWithStepInitFactory.h"
#include "Local_ODE_Xi_immob.h"


template <class T_DISCRETE_SYSTEM, class T_LINEAR_SOLVER>
class FunctionConcentrations
	: public ProcessLib::Process
{
public:
    enum In { Velocity=0 };
    enum Out { Concentrations = 0 };

	typedef FunctionConcentrations MyFunctionData; 
    typedef T_DISCRETE_SYSTEM MyDiscreteSystem;
    typedef T_LINEAR_SOLVER MyLinearSolver;

	// local matrix and vector
	typedef MathLib::LocalMatrix LocalMatrix;
	typedef MathLib::LocalVector LocalVector;

    // memory for discretized concentration vector
    typedef typename FemLib::FemNodalFunctionVector<MyDiscreteSystem>::type MyNodalFunctionVector;
	typedef typename FemLib::FemNodalFunctionScalar<MyDiscreteSystem>::type MyNodalFunctionScalar;

    // local assembler
	// for the linear systems, use the same settings as Mass_Transport
    typedef LinearTransportTimeODELocalAssembler<NumLib::ElementWiseTimeEulerEQSLocalAssembler> MyLinearAssemblerType;          
    typedef LinearTransportTimeODELocalAssembler<NumLib::ElementWiseTimeEulerResidualLocalAssembler> MyLinearResidualAssemblerType; 
    typedef LinearTransportJacobianLocalAssembler MyLinearJacobianAssemblerType;                                                      
	// for the nonlinear part, use different settings
	typedef NonLinearReactiveTransportTimeODELocalAssembler<NumLib::ElementWiseTimeEulerEQSLocalAssembler, MyNodalFunctionScalar>      MyNonLinearAssemblerType; 
	typedef NonLinearReactiveTransportTimeODELocalAssembler<NumLib::ElementWiseTimeEulerResidualLocalAssembler, MyNodalFunctionScalar> MyNonLinearResidualAssemblerType; 
	typedef NonLinearReactiveTransportJacobianLocalAssembler<MyNodalFunctionScalar, MyFunctionData>                                    MyNonLinearJacobianAssemblerType; 
	typedef NestedOdeNRIterationStepInitializer<MyNodalFunctionScalar, MyFunctionData>                                                 MyNRIterationStepInitializer;
	typedef NumLib::DiscreteNRSolverWithStepInitFactory<MyNRIterationStepInitializer>                                                  MyDiscreteNonlinearSolverFactory; 
	
	// linear Equation definition
    typedef SolutionLib::TemplateFemEquation<
            MyDiscreteSystem,
            MyLinearSolver,
            MyLinearAssemblerType,
            MyLinearResidualAssemblerType,
            MyLinearJacobianAssemblerType
            > MyLinearEquationType;
    // linear IVBV problem definition
    typedef SolutionLib::FemIVBVProblem<
            MyDiscreteSystem,
            MyLinearEquationType
            > MyLinearTransportProblemType;
    // linear equation solution algorithm definition
    typedef SolutionLib::SingleStepFEM<
            MyLinearTransportProblemType,
            MyLinearSolver
            > MyLinearSolutionType;

	// nonlinear coupled part solving xi
	// Equation definition
	typedef SolutionLib::TemplateFemEquation<
		    MyDiscreteSystem,
			MyLinearSolver,
			MyNonLinearAssemblerType,
			MyNonLinearResidualAssemblerType,
			MyNonLinearJacobianAssemblerType
	        > MyNonLinearEquationType;
	// IVBV problem definition
	typedef SolutionLib::FemIVBVProblem<
		    MyDiscreteSystem,
	        MyNonLinearEquationType
	        > MyNonLinearReactiveTransportProblemType;
	// Solution algorithm definition
	typedef SolutionLib::SingleStepFEM<
			MyNonLinearReactiveTransportProblemType,
			MyLinearSolver, 
			MyDiscreteNonlinearSolverFactory
			> MyNonLinearSolutionType;

	// the general reduction problem part
	typedef SolutionLib::FemKinReduction<MyDiscreteSystem> MyKinReductionProblemType; 
	typedef SolutionLib::SingleStepKinReduction<MyFunctionData, 
		                                        MyKinReductionProblemType,
		                                        MyLinearTransportProblemType, 
		                                        MyLinearSolutionType, 
												MyNonLinearReactiveTransportProblemType, 
	                                            MyNonLinearSolutionType> MyKinReductionSolution; 

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
        BaseLib::releaseObject(_ReductionKin);
		BaseLib::releaseObject(_local_ode_xi_immob); 

		size_t i; 
		for (i=0; i < _concentrations.size(); i++)
			BaseLib::releaseObject(_concentrations[i]);
		for (i=0; i < _eta_mob.size(); i++)
		{
			BaseLib::releaseObject(_linear_problems[i]); 
			BaseLib::releaseObject(_linear_solutions[i]); 
	        BaseLib::releaseObject(_eta_mob[i]); 
		}
	    for (i=0; i < _eta_immob.size(); i++)
	        BaseLib::releaseObject(_eta_immob[i]);
		for (i=0; i < _xi_mob.size(); i++)
			BaseLib::releaseObject(_xi_mob[i]); 
        for (i=0; i < _xi_immob.size(); i++)
		{
			BaseLib::releaseObject(_xi_immob[i]); 
			BaseLib::releaseObject(_xi_immob_new[i]); 
		}

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
		INFO("Solving %s...", getProcessName().c_str());
        initializeTimeStep(time);
        getSolution()->solveTimeStep(time);
        updateOutputParameter(time);
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
        getSolution()->accept(time);
    };

	// set values function
	void set_eta_mob_node_values  ( size_t eta_mob_idx,   MyNodalFunctionScalar* new_eta_mob_node_values   ); 
	void set_eta_immob_node_values( size_t eta_immob_idx, MyNodalFunctionScalar* new_eta_immob_node_values ); 
	void set_xi_mob_node_values   ( size_t xi_mob_idx,    MyNodalFunctionScalar* new_xi_mob_node_values    ); 
	void set_xi_immob_node_values ( size_t xi_immob_idx,  MyNodalFunctionScalar* new_xi_immob_node_values  ); 

	// calculate the rates on each node
	void update_node_kin_reaction_rates(void); 

	// update the change of kinetic reaction rate over the change of xi_mob
	void update_node_kin_reaction_drates_dxi(void); 

	// calculate nodal local ode problem of xi_immob
	void calc_nodal_xi_immob_ode(double dt);
	
protected:
    virtual void initializeTimeStep(const NumLib::TimeStep &time);

    virtual void updateOutputParameter(const NumLib::TimeStep &time);

    virtual MyKinReductionSolution* getSolution() {return _solution;};

    virtual void output(const NumLib::TimeStep &time);

private:
    DISALLOW_COPY_AND_ASSIGN(FunctionConcentrations);

private:
	virtual void convert_conc_to_eta_xi(void); 
	virtual void convert_eta_xi_to_conc(void); 

    // linear problem and solution pointer
	std::vector<MyLinearTransportProblemType*> _linear_problems;
    std::vector<MyLinearSolutionType*>         _linear_solutions;

	// nonlinear equation, problem and solution pointer
	MyNonLinearEquationType* _non_linear_eqs; 
	MyNonLinearReactiveTransportProblemType* _non_linear_problem;
	MyNonLinearSolutionType* _non_linear_solution; 
	Local_ODE_Xi_immob* _local_ode_xi_immob; 
	
	// reduction problem and solution
	MyKinReductionProblemType* _problem; 
	MyKinReductionSolution*    _solution; 

    FemLib::LagrangianFeObjectContainer* _feObjects; 
    
	NumLib::DiscreteDataConvergenceCheck _checker;

    // pointer to the reduction scheme. 
	ogsChem::chemReductionKin* _ReductionKin; 

    // concentrations vector, including all components in the MCP data structure
	std::vector<MyNodalFunctionScalar*> _concentrations; 
    // eta vector, including eta_mobile and eta_immobile
    std::vector<MyNodalFunctionScalar*> _eta_mob; 
	std::vector<MyNodalFunctionScalar*> _eta_immob;
    // xi_mobile and xi_immobile
    std::vector<MyNodalFunctionScalar*> _xi_mob; 
    std::vector<MyNodalFunctionScalar*> _xi_immob; 
	std::vector<MyNodalFunctionScalar*> _xi_immob_new; 
	// nodal reaction rates
	std::vector<MyNodalFunctionScalar*> _xi_mob_rates;
	std::vector<MyNodalFunctionScalar*> _xi_mob_drates_dxi; 
    std::vector<MyNodalFunctionScalar*> _xi_immob_rates;
}; 

#include "Concentrations.hpp"

#endif  // end of ifndef