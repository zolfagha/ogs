/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Concentration.h
 *
 * Created on 24.05.2013 by Reza Zolfaghari and Haibing Shao
 */

#include <functional>
#include "logog.hpp"

#include "MathLib/DataType.h"
#include "DiscreteLib/Utils/DiscreteSystemContainerPerMesh.h"
#include "FemLib/Function/FemNodalFunction.h"
#include "NumLib/Function/TXFunctionBuilder.h"
#include "NumLib/Function/TXFunctionDirect.h"
#include "OutputIO/OutputBuilder.h"
#include "OutputIO/OutputTimingBuilder.h"
#include "SolutionLib/Fem/FemSourceTerm.h"
#include "Ogs6FemData.h"
#include "NestedGIALocalProbNRIterationStepInitializer.h"
#include "MathLib/ODE/RungeKutta4.h"
#include "StepperBulischStoer.h"
#include "NonLinearGIATimeODELocalAssembler.h"

template <class T1, class T2>
bool FunctionReductConc<T1,T2>::initialize(const BaseLib::Options & option)
{
	size_t i;  // index
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    _msh_id = option.getOptionAsNum<size_t>("MeshID");
    size_t time_id = option.getOptionAsNum<size_t>("TimeGroupID");
    _tim = femData->list_tim[time_id];

    //mesh and FE objects
    MeshLib::IMesh* msh = femData->list_mesh[_msh_id];
    MyDiscreteSystem* dis = 0;
    dis = DiscreteLib::DiscreteSystemContainerPerMesh::getInstance()->createObject<MyDiscreteSystem>(msh);
    _feObjects = new FemLib::LagrangeFeObjectContainer(msh);

	// get the transformation class instance here
	this->_ReductionGIA = femData->m_GIA_ReductScheme;
    // make sure the reduction scheme is already initialized.
	if ( !(this->_ReductionGIA->IsInitialized()) )
	{
		// error msg
	    ERR("While initialize the Global Implicit Reactive Transport Process, the reduction class has not been correctly initialized! ");
		// then stop the program
		exit(1);
	}

	// first get the number of components
	_n_Comp = femData->map_ChemComp.size();
    // also get the size of secondary variables
    _n_eta   				= _ReductionGIA->get_n_eta  ();
    _n_eta_bar 				= _ReductionGIA->get_n_eta_bar();
    _n_xi_global    		= _ReductionGIA->get_n_xi_global ();
    _n_xi_local 			= _ReductionGIA->get_n_xi_local ();
	_n_xi_Sorp_tilde        = _ReductionGIA->get_n_xi_Sorp_tilde();
	_n_xi_Min_tilde         = _ReductionGIA->get_n_xi_Min_tilde();
	_n_xi_Sorp				= _ReductionGIA->get_n_xi_Sorp();
	_n_xi_Min				= _ReductionGIA->get_n_xi_Min();
	_I_mob					= _ReductionGIA->get_n_Comp_mob();
	_I_sorp                 = _ReductionGIA->get_n_Comp_sorb(); 
	_I_min					= _ReductionGIA->get_n_Comp_min();
	_I_kin                  = _ReductionGIA->get_n_Comp_kin();
	_n_xi_Kin_bar           = _ReductionGIA->get_n_xi_Kin_bar();
	_n_xi_Mob				= _ReductionGIA->get_n_xi_Mob();
	_n_xi_Sorp_bar			= _ReductionGIA->get_n_xi_Sorp_bar();
	_n_xi_Min_bar			= _ReductionGIA->get_n_xi_Min_bar();
    _n_xi_Kin               = _ReductionGIA->get_n_xi_Kin();
    _n_xi_Kin_bar           = _ReductionGIA->get_n_xi_Kin_bar();
    _J_tot_kin 				= _ReductionGIA->get_n_xi_Kin_total();

	// set concentrations of all components as output
	for ( i=0; i < _n_Comp; i++ )
		this->setOutputParameterName( i, femData->map_ChemComp[i]->second->get_name() );

	// creating local memory space to store IC and BC
	// initialize eta_mob,
	for ( i=0; i < _n_eta ; i++)
	{
		MyNodalFunctionScalar* eta_i = new MyNodalFunctionScalar();
		eta_i->initialize( *dis, FemLib::PolynomialOrder::Linear, 0.0 );
	    _eta.push_back(eta_i);
	}
	// initialize eta_immob,
	for ( i=0; i < _n_eta_bar; i++)
	{
		MyNodalFunctionScalar* eta_i = new MyNodalFunctionScalar();
	    eta_i->initialize( *dis, FemLib::PolynomialOrder::Linear, 0.0 );
		_eta_bar.push_back(eta_i);
	}
	// initialize xi_global
	for ( i=0; i < _n_xi_global ; i++ )
	{
		MyNodalFunctionScalar* xi_global_tmp1 = new MyNodalFunctionScalar();  
        MyNodalFunctionScalar* xi_global_tmp2 = new MyNodalFunctionScalar();  
		xi_global_tmp1->initialize( *dis, FemLib::PolynomialOrder::Linear, 0.0  );
        xi_global_tmp2->initialize( *dis, FemLib::PolynomialOrder::Linear, 0.0  );
		_xi_global_pre.push_back(xi_global_tmp1);
        _xi_global_cur.push_back(xi_global_tmp2);
	}

	// initialize xi_local
	for ( i=0; i < _n_xi_local ; i++ )
	{
		MyNodalFunctionScalar* xi_local_tmp       = new MyNodalFunctionScalar();  // xi_immob
		MyNodalFunctionScalar* xi_local_new_tmp   = new MyNodalFunctionScalar();  // xi_mob_rates
		xi_local_tmp->initialize( *dis, FemLib::PolynomialOrder::Linear, 0.0  );
		xi_local_new_tmp->initialize(*dis, FemLib::PolynomialOrder::Linear, 0.0  );
		_xi_local_old.push_back(xi_local_tmp);
		_xi_local_new.push_back(xi_local_new_tmp);
	}

	//initialize _global_vec_Rate
	for (i= 0; i < _J_tot_kin; i++)
	{
		MyNodalFunctionScalar* global_vec_Rate_tmp = new MyNodalFunctionScalar();  // xi_global
        global_vec_Rate_tmp->initialize( *dis, FemLib::PolynomialOrder::Linear, 0.0  );
     	_global_vec_Rate.push_back(global_vec_Rate_tmp);
	}
	// initialize the rates for xi_sorp, xi_min and xi_kin
	for (i = 0; i < _n_xi_Sorp; i++)
	{
		MyNodalFunctionScalar* xi_sorp_rates = new MyNodalFunctionScalar();
		xi_sorp_rates->initialize(*dis, FemLib::PolynomialOrder::Linear, 0.0);
		_xi_sorp_rates.push_back(xi_sorp_rates);
	}
	for (i = 0; i < _n_xi_Min; i++)
	{
		MyNodalFunctionScalar* xi_min_rates = new MyNodalFunctionScalar();
		xi_min_rates->initialize(*dis, FemLib::PolynomialOrder::Linear, 0.0);
		_xi_min_rates.push_back(xi_min_rates);
	}
	for (i = 0; i < _n_xi_Kin; i++)
	{
		MyNodalFunctionScalar* xi_kin_rates = new MyNodalFunctionScalar();
		xi_kin_rates->initialize(*dis, FemLib::PolynomialOrder::Linear, 0.0);
		_xi_kin_rates.push_back(xi_kin_rates);
	}

	// initialize drates_dxi
	for ( i=0; i < _J_tot_kin * (_n_xi_local + _n_xi_global) ; i++ )
	{
		MyNodalFunctionScalar* drate_dxi_tmp = new MyNodalFunctionScalar();  // drate_dxi instances
		drate_dxi_tmp->initialize(       *dis, FemLib::PolynomialOrder::Linear, 0.0  );
		_drates_dxi.push_back(drate_dxi_tmp);
	}

	// linear assemblers
    MyLinearAssemblerType* linear_assembler = new MyLinearAssemblerType(_feObjects);
    MyLinearResidualAssemblerType* linear_r_assembler = new MyLinearResidualAssemblerType(_feObjects);
    MyLinearJacobianAssemblerType* linear_j_eqs = new MyLinearJacobianAssemblerType(_feObjects);

	// for the linear transport problem, variables are eta_mobile
	for ( i=0; i < _n_eta ; i++ )
	{
		// set up problem
		MyLinearTransportProblemType* linear_problem = new MyLinearTransportProblemType(dis);
		MyLinearEquationType* linear_eqs = linear_problem->createEquation();
		linear_eqs->initialize(linear_assembler, linear_r_assembler, linear_j_eqs);
		linear_problem->setTimeSteppingFunction(*_tim);

		// set up variables
		// in this case, the variables includes: eta_0, eta_1, eta_2......,
		// which are the concentrations of all eta_mobile
		std::stringstream str_tmp;
		str_tmp << "eta_" << i ;
		linear_problem->addVariable( str_tmp.str() );
		_linear_problems.push_back(linear_problem);
	}

	// define nonlinear problem
	_non_linear_problem = new MyNonLinearReactiveTransportProblemType(dis);
	_non_linear_eqs     = _non_linear_problem->createEquation();
	//_non_linear_eqs->initialize( non_linear_j_assembler );
	_non_linear_problem->setTimeSteppingFunction(*_tim);
	// for nonlinear coupled transport problem, variables are xi_mobile species
	for ( i=0; i < _n_xi_global ; i++ )
	{
		std::stringstream str_tmp;
		str_tmp << "xi_global_" << i ;
		_non_linear_problem->addVariable( str_tmp.str() );
	}

	// reduction problem
	_problem = new MyGIAReductionProblemType( dis, _ReductionGIA );
	_problem->setTimeSteppingFunction(*_tim);  // applying the same time stepping function for all linear, non-linear and reduction problems
	// add variables to the KinReduction problem class
	// for the KinReduction problem, variables are the concentrations of all chemical components
	// add all concentrations to discretized memory space
	for ( i=0; i < _n_Comp; i++ )
	{
		MyVariableConc* comp_conc = _problem->addVariable( femData->map_ChemComp[i]->second->get_name() );
        FemVariableBuilder var_builder;
        var_builder.doit(femData->map_ChemComp[i]->second->get_name(), option, msh, femData->geo, femData->geo_unique_name, _feObjects, comp_conc);
	}

	// backward flushing global vector of concentrations
	MyNodalFunctionScalar* tmp_conc;
	for (i=0; i < _n_Comp; i++)
	{
		SolutionLib::FemIC* femIC = _problem->getVariable(i)->getIC();
	    tmp_conc = new MyNodalFunctionScalar();
		if ( femIC )
		{
			// FemIC vector is not empty
			tmp_conc->initialize(*dis, _problem->getVariable(i)->getCurrentOrder(), 0.0);
			femIC->setup(*tmp_conc);
		}
		else
		{
			// FemIC vector is empty
			// initialize the vector with zeros
			tmp_conc->initialize(*dis, _problem->getVariable(i)->getCurrentOrder(), 1E-50);  //RZ: avoid zero. 
		}
		_concentrations.push_back( tmp_conc );
	}

	// convert IC _concentrations to eta and xi
	convert_conc_to_eta_xi();

	// set IC for eta mobile
	for ( i=0; i < _n_eta; i++ )
	{
		SolutionLib::FemIC* eta_ic = new SolutionLib::FemIC(msh);
		eta_ic->addDistribution( femData->geo->getDomainObj(), new NumLib::TXFunctionDirect<double>( _eta[i]->getDiscreteData() ) );
		_linear_problems[i]->getVariable(0)->setIC( eta_ic );  // be careful, here each linear problem only has one variable "eta".
	}
	// set IC for xi_global
	for ( i=0; i < _n_xi_global; i++ )
	{
		SolutionLib::FemIC* xi_global_ic = new SolutionLib::FemIC(msh);
		xi_global_ic->addDistribution( femData->geo->getDomainObj(), new NumLib::TXFunctionDirect<double>( _xi_global_pre[i]->getDiscreteData() ) );
		_non_linear_problem->getVariable(i)->setIC( xi_global_ic );
	}


    // set up linear solution
	for ( i=0; i < _n_eta; i++ )
	{
		MyLinearSolutionType* linear_solution = new MyLinearSolutionType( dis, this->_linear_problems[i] );
		MyLinearSolver* linear_solver = linear_solution->getLinearEquationSolver();
		const BaseLib::Options* optNum = option.getSubGroup("Numerics");
		linear_solver->setOption(*optNum);
		this->_linear_solutions.push_back( linear_solution );
	}

	// set up non-linear solution
	myNRIterator 	 = new MyNRIterationStepInitializer(*this);
	myNSolverFactory = new MyDiscreteNonlinearSolverFactory( myNRIterator );
	this->_non_linear_solution = new MyNonLinearSolutionType( dis, this->_non_linear_problem, this, myNSolverFactory );
    this->_non_linear_solution->getDofEquationIdTable()->setNumberingType(DiscreteLib::DofNumberingType::BY_POINT);  // global order
    this->_non_linear_solution->getDofEquationIdTable()->setLocalNumberingType(DiscreteLib::DofNumberingType::BY_VARIABLE);  // local order
	const BaseLib::Options* optNum = option.getSubGroup("Numerics");

    // linear solver
	MyLinearSolver* linear_solver = this->_non_linear_solution->getLinearEquationSolver();
	linear_solver->setOption(*optNum);
	// set nonlinear solver options
	this->_non_linear_solution->getNonlinearSolver()->setOption(*optNum);
    // get the nonlinear solution dof manager
    this->_nl_sol_dofManager = this->_non_linear_solution->getDofEquationIdTable();

	// set up solution
    _solution = new MyGIAReductionSolution(dis, _problem, this, _linear_problems, _linear_solutions, _non_linear_problem, _non_linear_solution);

    this->setOutput(Concentrations, _solution->getCurrentSolution(0));


#ifdef _DEBUG
    // -----------debugging, output eta and xi----------------------
    for (i=0; i < _n_eta; i++) {
        std::stringstream str_tmp;
		str_tmp << "eta_" << i ;
        OutputVariableInfo var1(str_tmp.str(), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _eta[i]);
        femData->outController.setOutput(var1.name, var1);
    }
    for (i=0; i < _n_eta_bar; i++) {
        std::stringstream str_tmp;
		str_tmp << "eta_bar_" << i ;
        OutputVariableInfo var1(str_tmp.str(), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _eta_bar[i]);
        femData->outController.setOutput(var1.name, var1);
    }
    for (i=0; i < _n_xi_global; i++) {
        std::stringstream str_tmp1, str_tmp2;
		str_tmp1 << "xi_global_" << i ;
        OutputVariableInfo var1(str_tmp1.str(), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _xi_global_cur[i]);
        femData->outController.setOutput(var1.name, var1);
    }
    for (i=0; i < _n_xi_local; i++) {
        std::stringstream str_tmp1, str_tmp2;
		str_tmp1 << "xi_local_" << i ;
        OutputVariableInfo var1(str_tmp1.str(), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _xi_local_new[i]);
        femData->outController.setOutput(var1.name, var1);
    }
//    // -----------------end of debugging-----------------------------
#endif

    linear_solver = NULL;
    optNum = NULL;

    return true;
}

template <class T1, class T2>
void FunctionReductConc<T1, T2>::initializeTimeStep(const NumLib::TimeStep & time)
{
	// Update current time step
	this->_current_time_step = &time;

	size_t i;
    const NumLib::ITXFunction *vel = this->getInput<NumLib::ITXFunction>(Velocity);

	// set velocity for linear problem
	for ( i=0; i < _linear_problems.size(); i++ ) {
		_linear_problems[i]->getEquation()->getLinearAssembler()->setVelocity(vel);
	}
	// set velocity for nonlinear problem as well
	_non_linear_solution->getResidualFunction()->setVelocity(vel);
    _non_linear_solution->getDxFunction()->setVelocity(vel);
}

template <class T1, class T2>
void FunctionReductConc<T1, T2>::updateOutputParameter(const NumLib::TimeStep &/*time*/)
{ 
   // convert eta and xi back to concentrations
    convert_eta_xi_to_conc();
}

template <class T1, class T2>
void FunctionReductConc<T1, T2>::output(const NumLib::TimeStep &/*time*/)
{
    size_t i;

    // update data for output
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    this->_msh_id = this->_problem->getDiscreteSystem()->getMesh()->getID();
	// set the new output
	for (i=0; i<_concentrations.size(); i++) {
		OutputVariableInfo var1(this->getOutputParameterName(i), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _concentrations[i]);
        femData->outController.setOutput(var1.name, var1);
    }

#ifdef _DEBUG
    // -----------debugging, output eta and xi----------------------
    for (i=0; i < _n_eta; i++) {
        std::stringstream str_tmp;
		str_tmp << "eta_" << i ;
        OutputVariableInfo var1(str_tmp.str(), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _eta[i]);
        femData->outController.setOutput(var1.name, var1);
    }
    for (i=0; i < _n_eta_bar; i++) {
        std::stringstream str_tmp;
		str_tmp << "eta_bar_" << i ;
        OutputVariableInfo var1(str_tmp.str(), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _eta_bar[i]);
        femData->outController.setOutput(var1.name, var1);
    }
    for (i=0; i < _n_xi_global; i++) {
        std::stringstream str_tmp1, str_tmp2;
		str_tmp1 << "xi_global_" << i ;
        OutputVariableInfo var1(str_tmp1.str(), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _xi_global_cur[i]);
        femData->outController.setOutput(var1.name, var1);
    }
    for (i=0; i < _n_xi_local; i++) {
        std::stringstream str_tmp1, str_tmp2;
		str_tmp1 << "xi_local_" << i ;
        OutputVariableInfo var1(str_tmp1.str(), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _xi_local_new[i]);
        femData->outController.setOutput(var1.name, var1);
    }
    // -----------end of debugging----------------------------------
#endif
}

template <class T1, class T2>
void FunctionReductConc<T1, T2>::convert_conc_to_eta_xi(void)
{
	size_t node_idx, i;
    size_t _n_eta = this->_ReductionGIA->get_n_eta();
    size_t _n_eta_bar = this->_ReductionGIA->get_n_eta_bar();
    size_t _n_xi_global = this->_ReductionGIA->get_n_xi_global();
    size_t _n_xi_local = this->_ReductionGIA->get_n_xi_local();
	// only when the reduction scheme is fully initialized
	if ( this->_ReductionGIA->IsInitialized() )
	{
		// local vectors
		LocalVector loc_eta;
		LocalVector loc_eta_bar;
		LocalVector loc_xi_global;
		LocalVector loc_xi_local;
		LocalVector loc_conc;
		// allocate the memory for local vectors
		loc_eta     	= LocalVector::Zero( _n_eta );
		loc_eta_bar   	= LocalVector::Zero( _n_eta_bar );
		loc_xi_global   = LocalVector::Zero( _n_xi_global );
		loc_xi_local    = LocalVector::Zero( _n_xi_local );
		loc_conc        = LocalVector::Zero( _n_Comp );

		// for each nodes,
		for (node_idx=_concentrations[0]->getDiscreteData()->getRangeBegin();
			 node_idx < _concentrations[0]->getDiscreteData()->getRangeEnd();
			 node_idx++ )
		{
			for (i=0; i < _n_Comp; i++)
			{
				// gather all the concentrations
				loc_conc[i] = _concentrations[i]->getValue(node_idx);
			}  // end of for i

			// pass them to the transform function in the reductionKin class
			// and that the loc_eta_mob, local_eta_immob and local_xi
			this->_ReductionGIA->Conc2EtaXi( loc_conc, loc_eta, loc_eta_bar, loc_xi_global, loc_xi_local );

			// put the local eta and xi into the global vector
			// fill in eta_mob
			for (i=0; i < _n_eta; i++)
				this->_eta[i]->setValue(node_idx, loc_eta[i]);
			// fill in eta_immob
			for (i=0; i < _n_eta_bar; i++)
				this->_eta_bar[i]->setValue(node_idx, loc_eta_bar[i]);
			// fill in xi_mob
			for (i=0; i < _n_xi_global; i++)
			{
				//this->_xi_global_pre[i]->setValue(node_idx, loc_xi_global[i]);
				this->_xi_global_cur[i]->setValue(node_idx, loc_xi_global[i]);
				this->_xi_global_pre[i]->setValue(node_idx, loc_xi_global[i]);
			}

			for (i=0; i < _n_xi_local; i++)
            {
				// also over-write xi_local_new
                this->_xi_local_new[i]->setValue(node_idx, loc_xi_local[i]);
                this->_xi_local_old[i]->setValue(node_idx, loc_xi_local[i]);
            }
		}  // end of for node_idx

	}  // end of if _ReductionGIA
}

template <class T1, class T2>
void FunctionReductConc<T1, T2>::convert_eta_xi_to_conc(void)
{
	size_t node_idx, i;
    size_t _n_eta = this->_ReductionGIA->get_n_eta();
    size_t _n_eta_bar = this->_ReductionGIA->get_n_eta_bar();

	// only when the reduction scheme is fully initialized
	if ( this->_ReductionGIA->IsInitialized() )
	{
		// local vectors
		LocalVector loc_eta;
		LocalVector loc_eta_bar;
		LocalVector loc_xi_global;
		LocalVector loc_xi_local;
		LocalVector loc_conc;
		// allocate the memory for local vectors
		loc_eta		     = LocalVector::Zero( _n_eta );
		loc_eta_bar   	 = LocalVector::Zero( _n_eta_bar );
		loc_xi_global    = LocalVector::Zero( _n_xi_global );
		loc_xi_local     = LocalVector::Zero( _n_xi_local );
		loc_conc         = LocalVector::Zero( _n_Comp );

		// for each nodes,
		for (node_idx = _concentrations[0]->getDiscreteData()->getRangeBegin();
			 node_idx < _concentrations[0]->getDiscreteData()->getRangeEnd()  ;
			 node_idx++ )
		{
			// put the local eta and xi into the global vector
			// fill in eta_mob
			for (i=0; i < _n_eta; i++)
				loc_eta[i] = this->_eta[i]->getValue(node_idx);
			// fill in eta_immob
			for (i=0; i < _n_eta_bar; i++)
				loc_eta_bar[i] = this->_eta_bar[i]->getValue(node_idx);
			// fill in xi
			// this->_xi->setNodalValues( &loc_xi, node_idx*n_xi, n_xi );
			for (i=0; i < _n_xi_global; i++)
				loc_xi_global[i] = this->_xi_global_cur[i]->getValue(node_idx); // take the current time step one
		    for (i=0; i < _n_xi_local; i++)
				// using the xi_immob_new values
                loc_xi_local[i] = this->_xi_local_new[i]->getValue(node_idx);

			// pass them to the transform function in the reductionKin class
			// and that the loc_eta_mob, local_eta_immob and local_xi
			this->_ReductionGIA->EtaXi2Conc(loc_eta, loc_eta_bar, loc_xi_global, loc_xi_local, loc_conc);

			for (i=0; i < _n_Comp; i++)
			{
				// gather all the concentrations
				_concentrations[i]->setValue(node_idx, loc_conc[i]);
			}  // end of for i
		}  // end of for node_idx
	}  // end of if _ReductionGIA
}

template <class T1, class T2>
void FunctionReductConc<T1, T2>::set_eta_node_values( size_t eta_idx, MyNodalFunctionScalar* new_eta_node_values )
{
    size_t node_idx;
    for (node_idx = _concentrations[0]->getDiscreteData()->getRangeBegin();
	     node_idx < _concentrations[0]->getDiscreteData()->getRangeEnd();
		 node_idx++ )
        this->_eta[eta_idx]->setValue( node_idx, new_eta_node_values->getValue( node_idx ) );
}

template <class T1, class T2>
void FunctionReductConc<T1, T2>::set_eta_bar_node_values( size_t eta_bar_idx, MyNodalFunctionScalar* new_eta_bar_node_values )
{
    size_t node_idx;
    for (node_idx = _concentrations[0]->getDiscreteData()->getRangeBegin();
	     node_idx < _concentrations[0]->getDiscreteData()->getRangeEnd();
		 node_idx++ )
        this->_eta_bar[eta_bar_idx]->setValue( node_idx, new_eta_bar_node_values->getValue( node_idx ) );
}

template <class T1, class T2>
template <class T_X>
void FunctionReductConc<T1, T2>::update_xi_global_cur_nodal_values( const T_X & x_new )
{
    size_t n_var;
    n_var = this->_xi_global_cur.size();

    // T_X xi_mob_new_sol;
    // distribute solution vector to local vector for each variable
    for (size_t i=0; i<n_var; i++) {
        DiscreteLib::setLocalVector( *_nl_sol_dofManager, i, _msh_id, x_new, *this->_xi_global_cur[i]->getDiscreteData() );
    }

}


template <class T1, class T2>
void FunctionReductConc<T1, T2>::calc_nodal_local_problem(double dt, const double iter_tol, const double rel_tol, const double max_iter)
{
	size_t i, node_idx;
	double t0 = 0.0;
	double theta_water_content (0.5);  //monod2d example


	//pointer to the local problem
   // LocalProblem* _local_ode_xi_immob_GIA;
    _local_ode_xi_immob_GIA = new Local_ODE_Xi_immob_GIA( _ReductionGIA);

    MathLib::LocalMatrix mat_A2kin = _ReductionGIA->get_matrix_A2kin();

	//get the transformation matrices
	MathLib::LocalMatrix _mat_c_mob_2_xi_mob     = _ReductionGIA->get_matrix_C2Xi();
	MathLib::LocalMatrix _mat_c_immob_2_xi_immob = _ReductionGIA->get_matrix_Cbar2XiBar();

	// initialize the local vector
	MathLib::LocalVector loc_eta;
	MathLib::LocalVector loc_etabar;
	MathLib::LocalVector loc_xi_global;
	MathLib::LocalVector loc_xi_local;
	MathLib::LocalVector loc_xi_local_old_dt;
	MathLib::LocalVector loc_conc;
	MathLib::LocalVector vec_unknowns;
	MathLib::LocalVector loc_XiSorpTilde, loc_XiMinTilde, loc_XiKin, loc_XiBarKin, loc_XiBarKin_old;
	MathLib::LocalVector vec_conc, vec_XiBarKin, loc_xi_mobile, vec_xi_mob, loc_xi_immobile, vec_XiSorpBar, vec_XiMinBar, loc_xi_local_new, vec_unknowns_new;
	MathLib::LocalVector vec_conc_Mob, vec_conc_Sorp, vec_conc_Min, vec_conc_Kin, vec_conc_updated;
	// HS: notice that only concentrations of mob and sorp components 
	// will be converted to ln scale. 
	MathLib::LocalVector ln_conc_Mob, ln_conc_Sorp; 

	// initialize the local vector
	loc_eta        				= LocalVector::Zero( _n_eta );
	loc_etabar          		= LocalVector::Zero( _n_eta_bar );
	loc_xi_global       		= LocalVector::Zero( _n_xi_global );
	loc_xi_local       			= LocalVector::Zero( _n_xi_local );
	loc_xi_local_old_dt			= LocalVector::Zero( _n_xi_local );
	loc_conc	       			= LocalVector::Zero( _n_Comp );
	loc_XiSorpTilde				= LocalVector::Zero( _n_xi_Sorp_tilde);
	loc_XiMinTilde				= LocalVector::Zero( _n_xi_Min_tilde);
	loc_XiKin					= LocalVector::Zero( _n_xi_Kin);
	loc_XiBarKin				= LocalVector::Zero( _n_xi_Kin_bar);
	loc_XiBarKin_old			= LocalVector::Zero( _n_xi_Kin_bar);
	vec_unknowns				= LocalVector::Zero(_n_Comp ); 

	vec_conc            		= LocalVector::Zero(_n_Comp);
	vec_conc_Mob                = LocalVector::Zero(_I_mob);
	vec_conc_Sorp               = LocalVector::Zero(_I_sorp);
	vec_conc_Min                = LocalVector::Zero(_I_min);
	vec_conc_Kin                 = LocalVector::Zero(_I_kin);
	vec_conc_updated            = LocalVector::Zero(_n_Comp);

	ln_conc_Mob                 = LocalVector::Zero(_I_mob);
	ln_conc_Sorp                = LocalVector::Zero(_I_sorp);

	loc_xi_mobile		= LocalVector::Zero(_n_xi_Mob + _n_xi_Sorp + _n_xi_Min + _n_xi_Kin);
	loc_xi_immobile     = LocalVector::Zero(_n_xi_Sorp_bar + _n_xi_Min_bar + _n_xi_Kin_bar);

	vec_xi_mob			= LocalVector::Zero(_n_xi_Mob);
	vec_XiSorpBar		= LocalVector::Zero(_n_xi_Sorp_bar);
	vec_XiMinBar		= LocalVector::Zero(_n_xi_Min_bar);
	vec_XiBarKin        = LocalVector::Zero(_n_xi_Kin_bar);
	
	loc_xi_local_new	= LocalVector::Zero(_n_xi_local);
	vec_unknowns_new	= LocalVector::Zero(_n_Comp + _n_xi_Kin_bar);
	
	ogsChem::LocalVector  vec_xi_kin_rate  = ogsChem::LocalVector::Zero(_J_tot_kin);

	// signal, solving local ODEs of xi_immob
	INFO("--Solving local problem for xi_local and concentrations:");


    //  initialize the ODE using the StepperBulischStoer class
    //MathLib::StepperBulischStoer<Local_ODE_Xi_immob_GIA>* _sbs
	_sbs   = new MathLib::StepperBulischStoer<Local_ODE_Xi_immob_GIA>(loc_XiBarKin,
        															vec_xi_kin_rate,
        															t0,
        															1.0e-12,
        															1.0e-12,
        															true);

	//pointer to the local problem
    //LocalProblem* pSolve;
	_pSolve = new LocalProblem( _ReductionGIA, _sbs, _local_ode_xi_immob_GIA);

	// loop over all the nodes
    for (node_idx = _concentrations[0]->getDiscreteData()->getRangeBegin();
	     node_idx < _concentrations[0]->getDiscreteData()->getRangeEnd();
		 node_idx++ )
	{

			// on each node, get the right start value
			// get the right set of eta and xi
			for (i=0; i < _n_eta; i++)
				loc_eta[i] = this->_eta[i]->getValue(node_idx);
			// fill in eta_immob
			for (i=0; i < _n_eta_bar; i++)
				loc_etabar[i] = this->_eta_bar[i]->getValue(node_idx);
			for (i=0; i < _n_xi_global; i++)
				loc_xi_global[i] = this->_xi_global_cur[i]->getValue(node_idx); // using the current time step value
			for (i=0; i < _n_xi_local; i++)
			{
				loc_xi_local[i] 	   = this->_xi_local_new[i]->getValue(node_idx);
				loc_xi_local_old_dt[i] = this->_xi_local_old[i]->getValue(node_idx);
			}

			// get concentration values based on current updated eta and xi values. 
			this->_ReductionGIA->EtaXi2Conc(loc_eta, loc_etabar, loc_xi_global, loc_xi_local, loc_conc); 

			// skip the boundary nodes
			if ( ! this->_solution->isBCNode(node_idx) )
			{

			// need to record both old and new values of xi_kin
		    // for ODE solver
			loc_XiBarKin     = loc_xi_local.tail(_n_xi_Kin_bar);
			loc_XiBarKin_old = loc_xi_local_old_dt.tail(_n_xi_Kin_bar);

			vec_conc_Mob     = loc_conc.head(_I_mob); 
			vec_conc_Sorp    = loc_conc.segment(_I_mob, _I_sorp); 
			vec_conc_Min     = loc_conc.segment(_I_mob + _I_sorp, _I_min);
			vec_conc_Kin     = loc_conc.tail(_I_kin); 
			ln_conc_Mob.setZero();
			ln_conc_Sorp.setZero(); 
		    // convert the ln mobile conc to mobile conc
			_pSolve->cal_ln_conc_vec(_I_mob,  vec_conc_Mob,  ln_conc_Mob);
			_pSolve->cal_ln_conc_vec(_I_sorp, vec_conc_Sorp, ln_conc_Sorp);
			// now writting into the unknown vector. 
			vec_unknowns.head(_I_mob) = ln_conc_Mob;  // ln scale
			vec_unknowns.segment(_I_mob, _I_sorp) = ln_conc_Sorp;  // ln scale
			vec_unknowns.segment(_I_mob + _I_sorp, _I_min) = vec_conc_Min;  // linear scale
			vec_unknowns.tail(_I_kin) = vec_conc_Kin;  // linear scale

			// solve the local problem
			_pSolve->solve_LocalProblem_Newton_LineSearch( node_idx, 
				                                           dt, 
														   iter_tol, 
														   rel_tol, 
														   max_iter, 
														   vec_unknowns, 
														   loc_eta, 
														   loc_etabar, 
														   loc_xi_local, 
														   loc_xi_global, 
														   loc_XiBarKin_old);

		    // convert the ln mobile conc to mobile conc
			ln_conc_Mob  = vec_unknowns.head(_I_mob);
			ln_conc_Sorp = vec_unknowns.segment(_I_mob, _I_sorp); 
			vec_conc_Min = vec_unknowns.segment(_I_mob + _I_sorp, _I_min); 
			vec_conc_Kin = vec_unknowns.tail(_I_kin); 
			_pSolve->cal_exp_conc_vec(_I_mob, ln_conc_Mob, vec_conc_Mob);
			_pSolve->cal_exp_conc_vec(_I_sorp, ln_conc_Sorp, vec_conc_Sorp);

			vec_conc.head(_I_mob) = vec_conc_Mob; 
			vec_conc.segment(_I_mob, _I_sorp) = vec_conc_Sorp; 
			vec_conc.segment(_I_mob + _I_sorp, _I_min) = vec_conc_Min; 
			vec_conc.tail(_I_kin) = vec_conc_Kin; 			

			loc_xi_mobile     = _mat_c_mob_2_xi_mob * vec_conc_Mob;
		    vec_xi_mob        =  loc_xi_mobile.head(_n_xi_Mob);

		    loc_xi_immobile   = _mat_c_immob_2_xi_immob * vec_conc.tail(_I_sorp + _I_min + _I_kin);
		    vec_XiSorpBar     = loc_xi_immobile.head(_n_xi_Sorp_bar);
			vec_XiMinBar      = loc_xi_immobile.segment(_n_xi_Sorp_bar, _n_xi_Min_bar); 
			// directly take the value from ODE solution. 
			vec_XiBarKin      = _pSolve->get_vec_XiBarKin();

			// update the xi_local value. 
		    loc_xi_local_new.head   (_n_xi_Mob) 								 = vec_xi_mob;
		    loc_xi_local_new.segment(_n_xi_Mob, _n_xi_Sorp_bar)				     = vec_XiSorpBar;
		    loc_xi_local_new.segment(_n_xi_Mob + _n_xi_Sorp_bar, _n_xi_Min_bar)  = vec_XiMinBar; 
		    loc_xi_local_new.tail   (_n_xi_Kin_bar)							     = vec_XiBarKin; 

			vec_conc_updated.head(_I_mob) = vec_conc_Mob;
			vec_conc_updated.segment(_I_mob, _I_sorp) = vec_conc_Sorp;
		    vec_conc_updated.segment(_I_mob + _I_sorp, _I_min) = vec_conc_Min;
			vec_conc_updated.tail(_I_kin) = vec_conc_Kin;

			// collect the xi_local_new
			for (i=0; i < _n_xi_local; i++)
				_xi_local_new[i]->setValue(node_idx, loc_xi_local_new[i]);
			// collect new concentration values
			for (i=0; i < _n_Comp; i++)
				_concentrations[i]->setValue(node_idx, vec_conc_updated[i]);
		} // end of if
        else
        {
            // if it is boundary nodes, then xi_local_new is equal to xi_local.
            for (i=0; i < _n_xi_local; i++)
            	_xi_local_new[i]->setValue(node_idx, loc_xi_local[i] );

            // if it is boundary nodes, then conc_glob_new is equal to conc_glob.
            for (i=0; i < _n_Comp; i++)
            	_concentrations[i]->setValue(node_idx, loc_conc[i] );

        }  // end of else if
	}  // end of for node_idx

    //delete _pSolve; done in the header file
   // delete _sbs; done in the header file

}

template <class T1, class T2>
void FunctionReductConc<T1, T2>::set_BC_conc_node_values(std::size_t node_idx, std::size_t i_var, double node_value)
{
        _concentrations[i_var]->setValue( node_idx, node_value);
}

template <class T1, class T2>
void FunctionReductConc<T1, T2>::copy_cur_xi_global_to_pre(void)
{
    size_t n_var, node_idx;
    n_var = this->_xi_global_cur.size();

    for (size_t i=0; i<n_var; i++) {
        // loop over all the nodes
        for (node_idx = _xi_global_cur[i]->getDiscreteData()->getRangeBegin();
             node_idx < _xi_global_cur[i]->getDiscreteData()->getRangeEnd();
             node_idx++ )
        {
            _xi_global_pre[i]->setValue( node_idx, _xi_global_cur[i]->getValue( node_idx ) ); 
        }
    }
}

template <class T1, class T2>
void FunctionReductConc<T1, T2>::copy_cur_xi_local_to_pre(void)
{
    size_t node_idx;

    for (size_t i=0; i<_n_xi_local; i++) {
        // loop over all the nodes
        for (node_idx = _xi_local_new[i]->getDiscreteData()->getRangeBegin();
             node_idx < _xi_local_new[i]->getDiscreteData()->getRangeEnd();
             node_idx++ )
        {
        	_xi_local_old[i]->setValue( node_idx, _xi_local_new[i]->getValue( node_idx ) );
        }
    }
}
