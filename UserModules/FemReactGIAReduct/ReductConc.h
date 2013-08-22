/**
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ReductConc.h
 *
 * Created on    2013-05-17 by Reza Zolfaghari & Haibing Shao
 */

#ifndef REDUCT_CONC_H
#define REDUCT_CONC_H

#include "BaseLib/CodingTools.h"
#include "NumLib/TransientAssembler/ElementWiseTimeEulerEQSLocalAssembler.h"
#include "NumLib/TransientAssembler/ElementWiseTimeEulerResidualLocalAssembler.h"
#include "NumLib/Function/DiscreteDataConvergenceCheck.h"
#include "SolutionLib/Fem/FemIVBVProblem.h"
#include "SolutionLib/Fem/TemplateFemEquation.h"
#include "SolutionLib/Fem/SingleStepFEM.h"
#include "ProcessLib/AbstractTransientProcess.h"
#include "MathLib/DataType.h"
#include "NumLib/Nonlinear/DiscreteNRSolverWithStepInitFactory.h"
#include "FemGIARedSolution.h"
#include "SingleStepGIAReduction.h"
#include "NestedGIALocalProbNRIterationStepInitializer.h"
#include "LinearTransportTimeODELocalAssemberML.h"
#include "UserModules/FemKinReactGIA/LinearTransportJacobianLocalAssembler.h"
#include "NonLinearGIATimeODELocalAssembler.h"
#include "NonLinearGIAJacabianLocalAssembler.h"
#include "LocalProblem.h"
#include  "ogs6THMC/ChemLib/chemReductionGIA.h"
#include "TemplateFemEquation_GIAReduct.h"
#include "SingleStepGIAReductionNonlinear.h"


template <class T_DISCRETE_SYSTEM, class T_LINEAR_SOLVER >
class FunctionReductConc
	: public ProcessLib::Process
{
public:
	// input variable is velocity
    enum In { Velocity=0 };
	// no output variable
    enum Out { Concentrations=0 };

	// local matrix and vector
	typedef MathLib::LocalMatrix LocalMatrix;
	typedef MathLib::LocalVector LocalVector;

	typedef FunctionReductConc     MyFunctionData;    // FEM problem class
    typedef T_DISCRETE_SYSTEM      MyDiscreteSystem;  // Discretization
    typedef T_LINEAR_SOLVER        MyLinearSolver;    // linear solver
	
    // memory for discretized concentration vector
    typedef typename FemLib::FemNodalFunctionVector<MyDiscreteSystem>::type MyNodalFunctionVector;
	typedef typename FemLib::FemNodalFunctionScalar<MyDiscreteSystem>::type MyNodalFunctionScalar;

    // local assembler
	// for the linear systems, use the same settings as Mass_Transport
    typedef LinearTransportTimeODELocalAssemblerML<NumLib::ElementWiseTimeEulerEQSLocalAssembler> MyLinearAssemblerType;
    typedef LinearTransportTimeODELocalAssemblerML<NumLib::ElementWiseTimeEulerResidualLocalAssembler> MyLinearResidualAssemblerType;
    typedef LinearTransportJacobianLocalAssembler MyLinearJacobianAssemblerType;                                        

	// to be changed
    // typedef NonLinearGIAJacobianLocalAssembler<MyNodalFunctionScalar, MyFunctionData>                                    MyNonLinearJacobianAssemblerType;
	typedef NestedGIALocalProbNRIterationStepInitializer<MyNodalFunctionScalar, MyFunctionData>                          MyNRIterationStepInitializer;
	typedef NumLib::DiscreteNRSolverWithStepInitFactory<MyNRIterationStepInitializer>                                    MyDiscreteNonlinearSolverFactory;
	
	/**
      * linear Equation definition
      */
    typedef SolutionLib::TemplateFemEquation<
            MyDiscreteSystem,
            MyLinearSolver,
            MyLinearAssemblerType,
            MyLinearResidualAssemblerType,
            MyLinearJacobianAssemblerType
            > MyLinearEquationType;

    /**
      * linear IVBV problem definition
      */
    typedef SolutionLib::FemIVBVProblem<
            MyDiscreteSystem,
            MyLinearEquationType
            > MyLinearTransportProblemType;

    /** 
      * linear equation solution algorithm definition
      */
    typedef SolutionLib::SingleStepFEM<
            MyLinearTransportProblemType,
            MyLinearSolver
            > MyLinearSolutionType;

	/**
      * nonlinear coupled part solving xi
      * Equation definition
      */
	typedef TemplateFemEquation_GIAReduct<
		    MyDiscreteSystem,
            MyFunctionData,
			MyLinearSolver
	        > MyNonLinearEquationType;
	/**
      * FEM IVBV problem definition
      */
	typedef SolutionLib::FemIVBVProblem<
		    MyDiscreteSystem,
	        MyNonLinearEquationType
	        > MyNonLinearReactiveTransportProblemType;
	
    /**
	  * Solution algorithm definition
	  */
	typedef SingleStepGIAReductionNonlinear<
			MyNonLinearReactiveTransportProblemType,
			MyFunctionData,
			MyLinearSolver, 
			MyDiscreteNonlinearSolverFactory
			> MyNonLinearSolutionType;

    /**
	  * the general reduction problem part
	  */
	typedef SolutionLib::FemGIARedSolution<MyDiscreteSystem> MyGIAReductionProblemType; 

	typedef SolutionLib::SingleStepGIAReduction<MyFunctionData, 
		                                        MyGIAReductionProblemType,
		                                        MyLinearTransportProblemType, 
		                                        MyLinearSolutionType, 
												MyNonLinearReactiveTransportProblemType, 
	                                            MyNonLinearSolutionType> MyGIAReductionSolution;
    typedef typename MyGIAReductionProblemType::MyVariable MyVariableConc;

    FunctionReductConc() 
        : Process("REACT_GIA_REDUCT", 1, 1),
          _feObjects(0)
    {
        // set default parameter name
		ProcessLib::Process::setInputParameterName(Velocity, "Velocity");
    };
	
    /**
     * destructor, reclaim the memory
     */
    virtual ~FunctionReductConc()
    {
    	BaseLib::releaseObject(myNRIterator);
    	BaseLib::releaseObject(myNSolverFactory);

    	BaseLib::releaseObject(_problem);
    	BaseLib::releaseObject(_solution);

    	BaseLib::releaseObject(_ReductionGIA);
    	//BaseLib::releaseObject(_local_ode_xi_immob);

    	BaseLib::releaseObject(_non_linear_eqs);
    	BaseLib::releaseObject(_non_linear_problem);
    	BaseLib::releaseObject(_non_linear_solution);

    	BaseLib::releaseObjectsInStdVector(_linear_problems);
    	BaseLib::releaseObjectsInStdVector(_linear_solutions);

    	BaseLib::releaseObjectsInStdVector(_concentrations);
    	BaseLib::releaseObjectsInStdVector(_eta);
    	BaseLib::releaseObjectsInStdVector(_eta_bar);
    	BaseLib::releaseObjectsInStdVector(_xi_global);
    	BaseLib::releaseObjectsInStdVector(_xi_local);
    	//        BaseLib::releaseObjectsInStdVector(_rates);
    	BaseLib::releaseObjectsInStdVector(_kin_rates);
    	//BaseLib::releaseObjectsInStdVector(_xi_global_drates_dxi);
    	//BaseLib::releaseObjectsInStdVector(_xi_local);
    	//BaseLib::releaseObjectsInStdVector(_xi_local_new);
    	//BaseLib::releaseObjectsInStdVector(_xi_local_rates);

    	BaseLib::releaseObject(_nl_sol_dofManager);
    	BaseLib::releaseObject(_feObjects);

    };

    /**
	  * initialization of the problem class
	  */
    virtual bool initialize(const BaseLib::Options &option);

    /**
	  * finalize but nothing to do here
	  */
    virtual void finalize() {};

	/**
	  * returns the convergence checker
	  */
    virtual NumLib::IConvergenceCheck* getConvergenceChecker()
    {
        return &_checker;
    }

	/**
	  * function to solve the current time step
	  */
	virtual int solveTimeStep(const NumLib::TimeStep &time)
    {
		INFO("Solving %s...", getProcessName().c_str());
        initializeTimeStep(time);
        getSolution()->solveTimeStep(time);
        updateOutputParameter(time);
        return 0;
    }

    /**
	  * function to suggest the next time step
	  */
    virtual double suggestNext(const NumLib::TimeStep &time_current) 
    {
        return getSolution()->suggestNext(time_current); 
    }

    /**
	  * called when this problem is awake
	  */
    virtual bool isAwake(const NumLib::TimeStep &time) 
    { 
        return getSolution()->isAwake(time);
    }

    /**
	  * called when this time step is accepted
	  */
    virtual void accept(const NumLib::TimeStep &time)
    {
         output(time);
		 getSolution()->accept(time);
    };

    /**
	  * set function for eta and xi
	  */
	void set_eta_node_values     	( std::size_t eta_idx,   MyNodalFunctionScalar* new_eta_node_values   );
	void set_eta_bar_node_values    ( std::size_t eta_bar_idx, MyNodalFunctionScalar* new_eta_bar_node_values );
	void set_xi_global_node_values  ( std::size_t xi_global_idx,    MyNodalFunctionScalar* new_xi_global_node_values    );
	//void set_xi_local_node_values   ( std::size_t xi_local_idx,    MyNodalFunctionScalar* new_xi_local_node_values    );
    
    template <class T_X>
    void update_xi_global_nodal_values  ( const T_X & x_new );
	
    void update_xi_local_node_values ( void );

	/**
      * calculate the reaction rates on each node
      */
	void update_node_kin_reaction_rates(void);

	/**
      * update the change of kinetic reaction rate over the change of xi_mob
      */
	void update_node_kin_reaction_drates_dxi(void);

	/**
      * calculate nodal local problem
      */
	void calc_nodal_local_problem(double dt, const double iter_tol, const double rel_tol, const double max_iter);

	/**
      * calculate node based residual of global problem
      */
	void GlobalResidualAssembler(const NumLib::TimeStep & delta_t, const SolutionLib::SolutionVector & u_cur_xiglob, SolutionLib::SolutionVector & residual_global);

	//void assembly(const NumLib::TimeStep & delta_t, const std::size_t _n_xi, MathLib::LocalVector &global_u1_tilde, MathLib::LocalVector &global_u0_tilde, MathLib::LocalVector &global_u1,
	//		MathLib::LocalVector &global_u0, std::vector<MyNodalFunctionScalar*> _global_vec_LHS, std::vector<MyNodalFunctionScalar*> _global_vec_RHS);

	void assembly(	const NumLib::TimeStep & delta_t, const std::size_t _n_xi, const SolutionLib::SolutionVector & u_cur_xiglob,
					SolutionLib::SolutionVector & residual_global, std::size_t switchOn);

	void GlobalJacobianAssembler(const NumLib::TimeStep & delta_t, const SolutionLib::SolutionVector & u_cur_xiglob, SolutionLib::SolutionVector & Jacobian_global);

	void NumDiff(std::size_t & col,
			 	 ogsChem::LocalVector & delta_xi,
			 	 ogsChem::LocalVector & f,
			 	ogsChem::LocalVector & f_old,
			 	 ogsChem::LocalVector & unknown,
			 	 ogsChem::LocalVector & DrateDxi);
	void Vprime (MathLib::LocalVector & vec_conc,
				  MathLib::LocalMatrix & mat_vprime);

	ogsChem::chemReductionGIA* getReductionGIA(){return _ReductionGIA;}

	std::vector<MyNodalFunctionScalar*> & get_xi_global() {return _xi_global;}
	std::vector<MyNodalFunctionScalar*> & get_xi_local() {return _xi_local;}
	std::vector<MyNodalFunctionScalar*> & get_eta() {return _eta;}
	std::vector<MyNodalFunctionScalar*> & get_eta_bar() {return _eta_bar;}
	std::vector<MyNodalFunctionScalar*> & get_global_vec_Rate() {return _global_vec_Rate; }
	FemLib::LagrangeFeObjectContainer* get_feObjects(){return _feObjects;}
	NumLib::ITimeStepFunction* get_time_step() {return _tim;}

protected:
    virtual void initializeTimeStep(const NumLib::TimeStep &time);

    /**
      * this function is called to exchange output parameters
      */ 
    virtual void updateOutputParameter(const NumLib::TimeStep &time);

    /**
      * get the pointer of solution class for current problem
      */ 
    virtual MyGIAReductionSolution* getSolution() {return _solution;};

    /**
      * output the result of current solution
      */ 
    virtual void output(const NumLib::TimeStep &time);

private:
    DISALLOW_COPY_AND_ASSIGN(FunctionReductConc);

private:
    /**
      * convert nodal concentration values to eta and xi
      */ 
	virtual void convert_conc_to_eta_xi(void); 

    /**
      * convert nodal eta and xi values to concentrations
      */ 
	virtual void convert_eta_xi_to_conc(void); 

	DiscreteLib::DofEquationIdTable* _dofManager;

	MyDiscreteSystem *_dis_sys;

	NumLib::ITimeStepFunction* _tim;

	/**
      * linear problems
      */
	std::vector<MyLinearTransportProblemType*> _linear_problems;

    /**
      * linear solutions
      */
    std::vector<MyLinearSolutionType*>         _linear_solutions;

    /**
      * Nonlinear iterator
      */
	MyNRIterationStepInitializer*              myNRIterator;
	
    /**
      * Nonlinear solver factory
      */
    MyDiscreteNonlinearSolverFactory*          myNSolverFactory;
    
    /**
      * nonlinear equation
      */ 
	MyNonLinearEquationType*                  _non_linear_eqs;
    
    /**
      * non-linear problem
      */ 
	MyNonLinearReactiveTransportProblemType * _non_linear_problem;

    /**
      * non-linear solution
      */
	MyNonLinearSolutionType*                  _non_linear_solution;

	/**
      * the nested local ODE problem
      */ 
	//Local_ODE_Xi_immob*                       _local_ode_xi_immob; 
	
	/**
      * reduction problem
      */ 
	MyGIAReductionProblemType* _problem;

    /** 
      * reduction solution
      */
	MyGIAReductionSolution*    _solution;

	/**
      * FEM object
      */
    FemLib::LagrangeFeObjectContainer* _feObjects;
    
    double _theta;

    /**
      * velocity function
      */
    NumLib::ITXFunction* _vel;

	/**
      * convergence checker
      */
	NumLib::DiscreteDataConvergenceCheck _checker;

    /**
      * pointer to the reduction scheme. 
      */
	ogsChem::chemReductionGIA* _ReductionGIA;

    /**
      * concentrations vector
      * including all components in the MCP data structure
      */
	std::vector<MyNodalFunctionScalar*> _concentrations;


    /**
      * residual vector for the first two global equation 3.43 and  3.44
      */
	//std::vector<MyNodalFunctionScalar*> _residual_global_res43;
	//std::vector<MyNodalFunctionScalar*> _residual_global_res44;


    /**
      * nodal eta_mobile values
      */ 
    std::vector<MyNodalFunctionScalar*> _eta;

    /**
      * nodal eta_immobile values
      */ 
	std::vector<MyNodalFunctionScalar*> _eta_bar;
    
    /**
      * nodal xi_global values
      */ 
    std::vector<MyNodalFunctionScalar*> _xi_global;

    /**
      * nodal xi_local values
      */ 
    std::vector<MyNodalFunctionScalar*> _xi_local;

    /**
      * nodal concentration values
      */
    std::vector<MyNodalFunctionScalar*>  _conc_glob;

    /**
      * new nodal concentration values
      */
    std::vector<MyNodalFunctionScalar*>  _conc_glob_new;

    /**
      * nodal xi_immobile values
      */ 
	std::vector<MyNodalFunctionScalar*> _xi_local_new;
	
    /**
      * nodal reaction rates
      */ 
	std::vector<MyNodalFunctionScalar*> _kin_rates;

    /**
      * derivative of nodal xi_mobile rates over xi_mobile
      */ 
	std::vector<MyNodalFunctionScalar*> _xi_global_drates_dxi;

    /**
      * nodal xi_global reaction rates
      */
    std::vector<MyNodalFunctionScalar*> _xi_global_rates;

    /**
      * nodal xi_local reaction rates
      */ 
    std::vector<MyNodalFunctionScalar*> _xi_local_rates;

    /**
      * global reaction rate vector
      */
    std::vector<MyNodalFunctionScalar*> _global_vec_Rate;

    //temp global vectors
//    std::vector<MyNodalFunctionScalar*> _global_cur_xi_Sorp_tilde;
//    std::vector<MyNodalFunctionScalar*> _global_cur_xi_Min_tilde;
//    std::vector<MyNodalFunctionScalar*> _global_cur_xi_Sorp;
//    std::vector<MyNodalFunctionScalar*> _global_cur_xi_Min;
//    std::vector<MyNodalFunctionScalar*> _global_cur_xi_Kin;
//
//    std::vector<MyNodalFunctionScalar*> _global_u1_tilde;
//    std::vector<MyNodalFunctionScalar*> _global_u0_tilde;
//    std::vector<MyNodalFunctionScalar*> _global_u1;
//    std::vector<MyNodalFunctionScalar*> _global_u0;
//    std::vector<MyNodalFunctionScalar*> _global_vec_LHS;
//    std::vector<MyNodalFunctionScalar*> _global_vec_RHS;
//
//    std::vector<MyNodalFunctionScalar*> _global_pre_xi_Sorp_tilde;
//    std::vector<MyNodalFunctionScalar*> _global_pre_xi_Min_tilde;
//    std::vector<MyNodalFunctionScalar*> _global_pre_xi_Sorp;
//    std::vector<MyNodalFunctionScalar*> _global_pre_xi_Min;
//    std::vector<MyNodalFunctionScalar*> _global_pre_xi_Kin;
//
//    std::vector<MyNodalFunctionScalar*> _residual_global_res43;
//    std::vector<MyNodalFunctionScalar*> _residual_global_res44;
//    std::vector<MyNodalFunctionScalar*> _residual_global_res45;
//    std::vector<MyNodalFunctionScalar*> _residual_global_res46;
//    std::vector<MyNodalFunctionScalar*> _global_vec_Rate_45;
//    std::vector<MyNodalFunctionScalar*> _global_vec_Rate_46;
//    std::vector<MyNodalFunctionScalar*> _residual_global;
//    std::vector<MyNodalFunctionScalar*> _global_vec_LHS;
//    std::vector<MyNodalFunctionScalar*> _global_vec_RHS;
//    std::vector<MyNodalFunctionScalar*> _global_vec_LHS_sorp;
//    std::vector<MyNodalFunctionScalar*> _global_vec_RHS_sorp;
    /**
      * degree of freedom equation ID talbe for the nonlinear problem
      */ 
    DiscreteLib::DofEquationIdTable * _nl_sol_dofManager;

    /**
      * the id of the msh applied in this process
      */ 
    std::size_t _msh_id;

	std::size_t _n_xi_local,  _n_xi_global, _n_eta, _n_eta_bar; //TODO uninitialized

	//TODO no need to make them member variables
	std::size_t _n_xi_Mob, _n_xi_Sorp_tilde,_n_xi_Sorp, _n_xi_Sorp_bar, _n_xi_Sorp_bar_li, _n_xi_Sorp_bar_ld
			,_n_xi_Min, _n_xi_Min_tilde, _n_xi_Min_bar, _n_xi_Kin, _n_xi_Kin_bar, _I_mob, _I_sorp, _I_min, _I_NMin_bar, _J_tot_kin;
    std::size_t _n_Comp;
}; 

#include "ReductConc.hpp"

#endif  // end of ifndef
