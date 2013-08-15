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
#include "NonLinearGIATimeODELocalAssembler.h"
//#include "NonLinearGIAJacobianLocalAssembler.h"
#define _DEBUG

template <class T1, class T2>
bool FunctionReductConc<T1,T2>::initialize(const BaseLib::Options &option)
{
	size_t i;  // index
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    _msh_id = option.getOptionAsNum<size_t>("MeshID");
    size_t time_id = option.getOptionAsNum<size_t>("TimeGroupID");
    NumLib::ITimeStepFunction* tim = femData->list_tim[time_id];

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
    _n_eta   		= this->_ReductionGIA->get_n_eta  ();
	_n_eta_bar 		= this->_ReductionGIA->get_n_eta_bar();
	_n_xi_global    = this->_ReductionGIA->get_n_xi_global ();
	_n_xi_local 	= this->_ReductionGIA->get_n_xi_local ();

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
		MyNodalFunctionScalar* xi_global_tmp       = new MyNodalFunctionScalar();  // xi_global
        MyNodalFunctionScalar* xi_rates_tmp = new MyNodalFunctionScalar();  // R_kin rates
		xi_global_tmp->initialize(       *dis, FemLib::PolynomialOrder::Linear, 0.0  );
        xi_rates_tmp->initialize( *dis, FemLib::PolynomialOrder::Linear, 0.0  );
		_xi_global.push_back(xi_global_tmp);
        _kin_rates.push_back(xi_rates_tmp);
	}

//	// initialize drates_dxi
//	for ( i=0; i < _n_xi_global; i++ )
//	{
//		MyNodalFunctionScalar* drate_dxi_tmp = new MyNodalFunctionScalar();  // drate_dxi instances
//		drate_dxi_tmp->initialize(       *dis, FemLib::PolynomialOrder::Linear, 0.0  );
//		_drates_dxi.push_back(drate_dxi_tmp);
//	}
	// initialize xi_local
	for ( i=0; i < _n_xi_local ; i++ )
	{
		MyNodalFunctionScalar* xi_local_tmp       = new MyNodalFunctionScalar();  // xi_immob
        MyNodalFunctionScalar* rates_tmp = new MyNodalFunctionScalar();  // xi_mob_rates
		MyNodalFunctionScalar* xi_local_new_tmp   = new MyNodalFunctionScalar();  // xi_mob_rates
		xi_local_tmp->initialize( *dis, FemLib::PolynomialOrder::Linear, 0.0  );
        rates_tmp->initialize( *dis, FemLib::PolynomialOrder::Linear, 0.0  );
		xi_local_new_tmp->initialize(*dis, FemLib::PolynomialOrder::Linear, 0.0  );
		_xi_local.push_back(xi_local_tmp);
        _kin_rates.push_back(rates_tmp);
		_xi_local_new.push_back(xi_local_new_tmp);
	}

	// linear assemblers
    MyLinearAssemblerType* linear_assembler = new MyLinearAssemblerType(_feObjects);
    MyLinearResidualAssemblerType* linear_r_assembler = new MyLinearResidualAssemblerType(_feObjects);
    MyLinearJacobianAssemblerType* linear_j_eqs = new MyLinearJacobianAssemblerType(_feObjects);

    //To Do: for xi global
	MyNonLinearAssemblerType* non_linear_assembler = new MyNonLinearAssemblerType(_feObjects, this->_ReductionGIA );
	MyNonLinearResidualAssemblerType* non_linear_r_assembler = new MyNonLinearResidualAssemblerType(_feObjects, this->_ReductionGIA );
	MyNonLinearJacobianAssemblerType* non_linear_j_assembler = new MyNonLinearJacobianAssemblerType(_feObjects, this->_ReductionGIA, this);

	//_local_ode_xi_immob = new Local_ODE_Xi_immob( this->_ReductionGIA );

	// for the linear transport problem, variables are eta_mobile
	for ( i=0; i < _n_eta ; i++ )
	{
		// set up problem
		MyLinearTransportProblemType* linear_problem = new MyLinearTransportProblemType(dis);
		MyLinearEquationType* linear_eqs = linear_problem->createEquation();
		linear_eqs->initialize(linear_assembler, linear_r_assembler, linear_j_eqs);
		linear_problem->setTimeSteppingFunction(*tim);

		// set up variables
		// in this case, the variables includes: eta_0, eta_1, eta_2......,
		// which are the concentrations of all eta_mobile
		std::stringstream str_tmp;
		str_tmp << "eta_" << i ;
		linear_problem->addVariable( str_tmp.str() );
		_linear_problems.push_back(linear_problem);
	}

	// TO DO : xi global
	// define nonlinear problem
	_non_linear_problem = new MyNonLinearReactiveTransportProblemType(dis);
	_non_linear_eqs     = _non_linear_problem->createEquation();
	_non_linear_eqs->initialize( non_linear_assembler, non_linear_r_assembler, non_linear_j_assembler );
	_non_linear_problem->setTimeSteppingFunction(*tim);
	// for nonlinear coupled transport problem, variables are xi_mobile species
	for ( i=0; i < _n_xi_global ; i++ )
	{
		std::stringstream str_tmp;
		str_tmp << "xi_global_" << i ;
		_non_linear_problem->addVariable( str_tmp.str() );
	}

	// reduction problem
	_problem = new MyGIAReductionProblemType( dis, _ReductionGIA );
	_problem->setTimeSteppingFunction(*tim);  // applying the same time stepping function for all linear, non-linear and reduction problems
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
			tmp_conc->initialize(*dis, _problem->getVariable(i)->getCurrentOrder(), 0.0);
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
		xi_global_ic->addDistribution( femData->geo->getDomainObj(), new NumLib::TXFunctionDirect<double>( _xi_global[i]->getDiscreteData() ) );
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

	//TODO : xi global
	// set up non-linear solution
	// NRIterationStepInitializer
	myNRIterator = new MyNRIterationStepInitializer(non_linear_r_assembler, non_linear_j_assembler);
	myNSolverFactory = new MyDiscreteNonlinearSolverFactory( myNRIterator );
	this->_non_linear_solution = new MyNonLinearSolutionType( dis, this->_non_linear_problem, myNSolverFactory );
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

//    // set initial output parameter
//	for (i=0; i<_concentrations.size(); i++) {
//		OutputVariableInfo var1(this->getOutputParameterName(i), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _concentrations[i]);
//        femData->outController.setOutput(var1.name, var1);
 //   }

#ifdef _DEBUG
    // -----------debugging, output eta and xi----------------------
    for (i=0; i < _n_eta; i++) {
        std::stringstream str_tmp;
		str_tmp << "eta_" << i ;
        OutputVariableInfo var1(str_tmp.str(), _msh_id,  OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _eta[i]);
        femData->outController.setOutput(var1.name, var1);
    }
//    for (i=0; i < _n_eta_immob; i++) {
//        std::stringstream str_tmp;
//		str_tmp << "eta_immob_" << i ;
//        OutputVariableInfo var1(str_tmp.str(), _msh_id,  OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _eta_immob[i]);
//        femData->outController.setOutput(var1.name, var1);
//    }
//    for (i=0; i < _n_xi_mob; i++) {
//        std::stringstream str_tmp1, str_tmp2;
//		str_tmp1 << "xi_mob_" << i ;
//        OutputVariableInfo var1(str_tmp1.str(), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _xi_mob[i]);
//        femData->outController.setOutput(var1.name, var1);
//        str_tmp2 << "xi_mob_rate_" << i;
//        OutputVariableInfo var2(str_tmp2.str(), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _xi_mob_rates[i]);
//        femData->outController.setOutput(var2.name, var2);
//    }
    for (i=0; i < _n_xi_local; i++) {
        std::stringstream str_tmp1, str_tmp2;
		str_tmp1 << "xi_local_" << i ;
        OutputVariableInfo var1(str_tmp1.str(), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _xi_local[i]);
        femData->outController.setOutput(var1.name, var1);
//        str_tmp2 << "xi_immob_rate_" << i;
//        OutputVariableInfo var2(str_tmp2.str(), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _xi_immob_rates[i]);
//        femData->outController.setOutput(var2.name, var2);
    }

    for (i=0; i < _n_Comp; i++) {
        std::stringstream str_tmp1, str_tmp2;
		str_tmp1 << "conc_" << i ;
        OutputVariableInfo var1(str_tmp1.str(), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _concentrations[i]);
        femData->outController.setOutput(var1.name, var1);
    }
//    // -----------------end of debugging-----------------------------
#endif

    linear_solver = NULL;
    optNum = NULL;

    return true;
}

template <class T1, class T2>
void FunctionReductConc<T1, T2>::initializeTimeStep(const NumLib::TimeStep &/*time*/)
{
	size_t i;
    const NumLib::ITXFunction *vel = this->getInput<NumLib::ITXFunction>(Velocity);

	// set velocity for linear problem
	for ( i=0; i < _linear_problems.size(); i++ ) {
		_linear_problems[i]->getEquation()->getLinearAssembler()->setVelocity(vel);
	//	_linear_problems[i]->getEquation()->getResidualAssembler()->setVelocity(vel);
	//	_linear_problems[i]->getEquation()->getJacobianAssembler()->setVelocity(vel);
	}
	// set velocity for nonlinear problem as well
	_non_linear_problem->getEquation()->getLinearAssembler()->setVelocity(vel);
    _non_linear_problem->getEquation()->getResidualAssembler()->setVelocity(vel);
	_non_linear_problem->getEquation()->getJacobianAssembler()->setVelocity(vel);
 //   // set xi_mob_rates for non-linear problem
	//_non_linear_problem->getEquation()->getLinearAssembler()  ->set_xi_mob_rates( &_xi_mob_rates ); 
	//_non_linear_problem->getEquation()->getResidualAssembler()->set_xi_mob_rates( &_xi_mob_rates ); 
	//_non_linear_problem->getEquation()->getJacobianAssembler()->set_xi_mob_rates( &_xi_mob_rates ); 
 //   // set xi_immob_rates for non-linear problem
	//_non_linear_problem->getEquation()->getLinearAssembler()  ->set_xi_immob_rates( &_xi_immob_rates ); 
	//_non_linear_problem->getEquation()->getResidualAssembler()->set_xi_immob_rates( &_xi_immob_rates ); 
	//_non_linear_problem->getEquation()->getJacobianAssembler()->set_xi_immob_rates( &_xi_immob_rates ); 
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
//    for (i=0; i < _n_eta_immob; i++) {
//        std::stringstream str_tmp;
//		str_tmp << "eta_immob_" << i ;
//        OutputVariableInfo var1(str_tmp.str(), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _eta_immob[i]);
//        femData->outController.setOutput(var1.name, var1);
//    }
//    for (i=0; i < _n_xi_mob; i++) {
//        std::stringstream str_tmp1, str_tmp2;
//		str_tmp1 << "xi_mob_" << i ;
//        OutputVariableInfo var1(str_tmp1.str(), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _xi_mob[i]);
//        femData->outController.setOutput(var1.name, var1);
//        str_tmp2 << "xi_mob_rate_" << i; 
//        OutputVariableInfo var2(str_tmp2.str(), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _xi_mob_rates[i]);
//        femData->outController.setOutput(var2.name, var2);
//    }
//    for (i=0; i < _n_xi_immob; i++) {
//        std::stringstream str_tmp1, str_tmp2;
//		str_tmp1 << "xi_immob_" << i ;
//        OutputVariableInfo var1(str_tmp1.str(), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _xi_immob[i]);
//        femData->outController.setOutput(var1.name, var1);
//        str_tmp2 << "xi_immob_rate_" << i; 
//        OutputVariableInfo var2(str_tmp2.str(), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _xi_immob_rates[i]);
//        femData->outController.setOutput(var2.name, var2);
//    }
//    // -----------end of debugging----------------------------------
#endif
}

template <class T1, class T2>
void FunctionReductConc<T1, T2>::convert_conc_to_eta_xi(void)
{
	size_t node_idx, i;

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
				this->_xi_global[i]->setValue(node_idx, loc_xi_global[i]);

			for (i=0; i < _n_xi_local; i++)
            {
				this->_xi_local[i]->setValue(node_idx, loc_xi_local[i]);
                // also over-write xi_immob_new
                this->_xi_local_new[i]->setValue(node_idx, loc_xi_local[i]);
            }
		}  // end of for node_idx

	}  // end of if _ReductionGIA
}

template <class T1, class T2>
void FunctionReductConc<T1, T2>::convert_eta_xi_to_conc(void)
{
	size_t node_idx, i;

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
				loc_xi_global[i] = this->_xi_global[i]->getValue(node_idx);
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
void FunctionReductConc<T1, T2>::set_xi_global_node_values( size_t xi_global_idx, MyNodalFunctionScalar* new_xi_global_node_values )
{
    size_t node_idx;
    for (node_idx = _concentrations[0]->getDiscreteData()->getRangeBegin();
	     node_idx < _concentrations[0]->getDiscreteData()->getRangeEnd();
		 node_idx++ )
        this->_xi_global[xi_global_idx]->setValue( node_idx, new_xi_global_node_values->getValue( node_idx ) );
}

template <class T1, class T2>
template <class T_X>
void FunctionReductConc<T1, T2>::update_xi_global_nodal_values( const T_X & x_new )
{
    size_t n_var;
    n_var = this->_xi_global.size();

    // T_X xi_mob_new_sol;
    // distribute solution vector to local vector for each variable
    for (size_t i=0; i<n_var; i++) {
        DiscreteLib::setLocalVector( *_nl_sol_dofManager, i, _msh_id, x_new, *this->_xi_global[i]->getDiscreteData() );
    }

}

template <class T1, class T2>
void FunctionReductConc<T1, T2>::update_xi_local_node_values( void )
{
    size_t i, node_idx;
    for ( i=0; i < _xi_local.size() ; i++ )
    for (node_idx = _concentrations[0]->getDiscreteData()->getRangeBegin();
	     node_idx < _concentrations[0]->getDiscreteData()->getRangeEnd();
		 node_idx++ )
    	this->_xi_local[i]->setValue( node_idx, _xi_local_new[i]->getValue( node_idx ) );
}

template <class T1, class T2>
void FunctionReductConc<T1, T2>::update_node_kin_reaction_rates(void)
{
	size_t node_idx, i;

	// only when the reduction scheme is fully initialized
	if ( this->_ReductionGIA->IsInitialized() )
	{
		// local vectors
		LocalVector loc_eta;
		LocalVector loc_eta_bar;
		LocalVector loc_xi_global;
		LocalVector loc_xi_local;
		LocalVector loc_conc;
        LocalVector loc_xi_global_rates;
        LocalVector loc_xi_local_rates;

		// allocate the memory for local vectors
		loc_eta        		= LocalVector::Zero( _n_eta );
		loc_eta_bar     	= LocalVector::Zero( _n_eta_bar );
		loc_xi_global       = LocalVector::Zero( _n_xi_global );
		loc_xi_local       	= LocalVector::Zero( _n_xi_local );
		loc_conc           	= LocalVector::Zero( _n_Comp );
        loc_xi_global_rates = LocalVector::Zero( _n_xi_global );
        loc_xi_local_rates 	= LocalVector::Zero( _n_xi_local );

		// for each nodes,
		for (node_idx = _concentrations[0]->getDiscreteData()->getRangeBegin();
			 node_idx < _concentrations[0]->getDiscreteData()->getRangeEnd()  ;
			 node_idx++ )
		{
            if ( ! this->_solution->isBCNode(node_idx) )
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
                    loc_xi_global[i] = this->_xi_global[i]->getValue(node_idx);
                for (i=0; i < _n_xi_local; i++)
                    // loc_xi_immob[i] = this->_xi_immob[i]->getValue(node_idx);
                    loc_xi_local[i] = this->_xi_local_new[i]->getValue(node_idx); // using new rates instead! // HS 26.11.2012

                // calculate rates;
                this->_ReductionGIA->Calc_Kin_Rate( loc_eta,
                                                   loc_eta_bar,
                                                   loc_xi_global,
                                                   loc_xi_local,
                                                   loc_xi_global_rates,
                                                   loc_xi_local_rates );

                // write the rates into the global nodal vector
                for (i=0; i < _n_xi_global; i++)
                    this->_xi_global_rates[i]->setValue( node_idx, loc_xi_global_rates[i] );

                for (i=0; i < _n_xi_local; i++)
                    this->_xi_local_rates[i]->setValue( node_idx, loc_xi_local_rates[i] );
            }  // end of if boundary node
            else
            {
                // if is boundary node, then set rates to zero.
                for (i=0; i < _n_xi_global; i++)
                    this->_xi_global_rates[i]->setValue( node_idx, 0.0 );

                for (i=0; i < _n_xi_local; i++)
                    this->_xi_local_rates[i]->setValue( node_idx, 0.0 );
            }

		}  // end of for node_idx

 //       _non_linear_problem->getEquation()->getResidualAssembler()->set_xi_global_rates(   & _xi_global_rates );
 //       _non_linear_problem->getEquation()->getResidualAssembler()->set_xi_local_rates( & _xi_local_rates );

        _non_linear_problem->getEquation()->getJacobianAssembler()->set_xi_mob_rates(   & _xi_global_rates );
        _non_linear_problem->getEquation()->getJacobianAssembler()->set_xi_immob_rates( & _xi_local_rates );
	}  // end of if _ReductionGIA
}


template <class T1, class T2>
void FunctionReductConc<T1, T2>::update_node_kin_reaction_drates_dxi(void)
{
	const double epsilon  = 1.0e-14;
    double drates_dxi_tmp = 0.0;

    size_t node_idx, i, j;

	LocalVector loc_eta;
	LocalVector loc_eta_bar;
	LocalVector loc_xi_global;
	LocalVector loc_xi_global_tmp;  // used for xi increment
	LocalVector loc_xi_local;
	LocalVector loc_xi_global_rates_base;
    LocalVector loc_xi_global_rates_new;
    LocalVector loc_xi_local_rates;

	// initialize the local vector
	loc_eta     		      = LocalVector::Zero( _n_eta );
	loc_eta_bar         	  = LocalVector::Zero( _n_eta_bar );
	loc_xi_global             = LocalVector::Zero( _n_xi_global );
	loc_xi_local          	  = LocalVector::Zero( _n_xi_local );
    loc_xi_global_rates_base  = LocalVector::Zero( _n_xi_global );
    loc_xi_global_rates_new   = LocalVector::Zero( _n_xi_global );
    loc_xi_local_rates    	  = LocalVector::Zero( _n_xi_local );

	// loop over all nodes
	for (node_idx = _concentrations[0]->getDiscreteData()->getRangeBegin();
	     node_idx < _concentrations[0]->getDiscreteData()->getRangeEnd()  ;
		 node_idx++)
	{
        // if not boundary nodes, calculates its drate/dxi.
        if ( ! this->_solution->isBCNode(node_idx) )
        {
            // read the local values
            for (i=0; i < _n_eta; i++)
                loc_eta[i] = this->_eta[i]->getValue(node_idx);
            // fill in eta_immob
            for (i=0; i < _n_eta_bar; i++)
                loc_eta_bar[i] = this->_eta_bar[i]->getValue(node_idx);
            for (i=0; i < _n_xi_global; i++)
                loc_xi_global[i] = this->_xi_global[i]->getValue(node_idx);
            for (i=0; i < _n_xi_local; i++)
                loc_xi_local[i] = this->_xi_local_new[i]->getValue(node_idx);  // HS, use the new immob values after ODE calculation

//            // calculate rates;
//            this->_ReductionGIA->Calc_Xi_Rate( loc_eta,
//                                               loc_eta_bar,
//                                               loc_xi_global,
//                                               loc_xi_local,
//                                               loc_xi_global_rates_base,
//                                               loc_xi_local_rates );

            // loop over each xi_mob,
            for (j=0; j < _n_xi_global; j++)
            {
                // get a clean copy of the origiinal xi_mob
                loc_xi_global_tmp = loc_xi_global;
                // give an increment to the xi value
                if ( abs( loc_xi_global(j) ) >= 1.0e-16 )  // loc_xi_mob(j) != 0.0
                    loc_xi_global_tmp(j) += epsilon * abs( loc_xi_global(j) );
                else
                    loc_xi_global_tmp(j) += epsilon;

                // calculate the new rate
                //            // calculate rates;
                //            this->_ReductionGIA->Calc_Xi_Rate( loc_eta,
                //                                               loc_eta_bar,
                //                                               loc_xi_global,
                //                                               loc_xi_local,
                //                                               loc_xi_global_tmp,
                //                                               loc_xi_local_rates );
                // loop over all xi_mob
                for (i=0; i < _n_xi_global; i++)
                {
                    // calculate derivative
                    if ( abs( loc_xi_global(j) ) >= 1.0e-16 )  // loc_xi_mob(j) != 0.0
                        drates_dxi_tmp = ( loc_xi_global_rates_new(i) - loc_xi_global_tmp(i) ) / epsilon /  abs( loc_xi_global(j) );
                    else
                        drates_dxi_tmp = ( loc_xi_global_rates_new(i) - loc_xi_global_tmp(i) ) /  epsilon ;

                     _xi_global_drates_dxi[i*_n_xi_global+j]->setValue( node_idx, drates_dxi_tmp );
                }  // end of for i
            }  // end of for j
        }  // end of if boundary node
        else
        {
            // if boundary node, then set drates/dxi to zeros
            for (i=0; i < _n_xi_global; i++)
                for (j=0; j < _n_xi_global; j++)
                    _xi_global_drates_dxi[i*_n_xi_global+j]->setValue( node_idx, 0.0 );
        }  // end of else
	}  // end of for node_idx

    // setting the rates
    _non_linear_problem->getEquation()->getJacobianAssembler()->set_xi_global_drate_dxi( &_xi_global_drates_dxi );
}


template <class T1, class T2>
void FunctionReductConc<T1, T2>::calc_nodal_local_problem(double dt, const double iter_tol, const double rel_tol, const double max_iter)
{
	size_t i, node_idx;

    //pointer to the local problem
    LocalProblem* pSolve;
    pSolve = new LocalProblem( _ReductionGIA);

	//get the transformation matrices
	_mat_c_mob_2_xi_mob     = _ReductionGIA->get_matrix_C2Xi();
	_mat_c_immob_2_xi_immob = _ReductionGIA->get_matrix_Cbar2XiBar();
	_n_xi_Sorp_tilde        = _ReductionGIA->get_n_xi_Sorp_tilde();
	_n_xi_Min_tilde         = _ReductionGIA->get_n_xi_Min_tilde();
	_n_xi_Sorp				= _ReductionGIA->get_n_xi_Sorp();
	_n_xi_Min				= _ReductionGIA->get_n_xi_Min();
	_I_mob					= _ReductionGIA->get_n_Comp_mob();
	_I_NMin_bar             = _ReductionGIA->get_n_Comp_NMin_bar();
	_I_min					= _ReductionGIA->get_n_Comp_min();
	_n_xi_Kin_bar           = _ReductionGIA->get_n_xi_Kin_bar();
	_n_xi_Mob				= _ReductionGIA->get_n_xi_Mob();
	_n_xi_Sorp_bar			= _ReductionGIA->get_n_xi_Sorp_bar();
	_n_xi_Min_bar			= _ReductionGIA->get_n_xi_Min_bar();
    _n_xi_Kin               = _ReductionGIA->get_n_xi_Kin(); 
    _n_xi_Kin_bar           = _ReductionGIA->get_n_xi_Kin_bar();

	// initialize the local vector
	MathLib::LocalVector loc_eta;
	MathLib::LocalVector loc_etabar;
	MathLib::LocalVector loc_xi_global;
	MathLib::LocalVector loc_xi_local;
	MathLib::LocalVector loc_xi_local_old_dt;
	MathLib::LocalVector loc_conc;
	MathLib::LocalVector vec_tot_mass_constrain;
	MathLib::LocalVector vec_unknowns;
	MathLib::LocalVector loc_XiSorpTilde, loc_XiMinTilde, loc_XiKin, loc_XiBarKin, loc_XiBarKin_old;
	MathLib::LocalVector vec_conc, vec_XiBarKin, loc_xi_mobile, vec_xi_mob, loc_xi_immobile, vec_XiSorpBar, vec_XiMinBar, loc_xi_local_new, vec_unknowns_new;

	// initialize the local vector
	loc_eta        				= LocalVector::Zero( _n_eta );
	loc_etabar          		= LocalVector::Zero( _n_eta_bar );
	loc_xi_global       		= LocalVector::Zero( _n_xi_global );
	loc_xi_local       			= LocalVector::Zero( _n_xi_local );
	loc_xi_local_old_dt			= LocalVector::Zero( _n_xi_local );
	loc_conc	       			= LocalVector::Zero( _n_Comp );
	vec_tot_mass_constrain	    = LocalVector::Zero( _n_eta + _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Kin + _n_eta_bar + _n_xi_Kin_bar );
//	loc_conc	       			= LocalVector::Zero( _n_Comp + _n_xi_Kin_bar);
	loc_XiSorpTilde				= LocalVector::Zero( _n_xi_Sorp_tilde);
	loc_XiMinTilde				= LocalVector::Zero( _n_xi_Min_tilde);
	loc_XiKin					= LocalVector::Zero( _n_xi_Kin);
	loc_XiBarKin				= LocalVector::Zero( _n_xi_Kin_bar);
	loc_XiBarKin_old			= LocalVector::Zero( _n_xi_Kin_bar);
	vec_unknowns				= LocalVector::Zero(_n_Comp + _n_xi_Kin_bar);

	vec_conc     		= LocalVector::Zero(_n_Comp);
	vec_XiBarKin		= LocalVector::Zero(_n_xi_Kin_bar);
	loc_xi_mobile		= LocalVector::Zero(_n_xi_Mob + _n_xi_Sorp + _n_xi_Min + _n_xi_Kin);
	vec_xi_mob			= LocalVector::Zero(_n_xi_Mob);
	loc_xi_immobile		= LocalVector::Zero(_n_xi_Sorp_bar + _n_xi_Min_bar + _n_xi_Kin_bar);
	vec_XiSorpBar		= LocalVector::Zero(_n_xi_Sorp_bar);
	vec_XiMinBar		= LocalVector::Zero(_n_xi_Min_bar);
	loc_xi_local_new	= LocalVector::Zero(_n_xi_local);
	vec_unknowns_new	= LocalVector::Zero(_n_Comp + _n_xi_Kin_bar);


	// signal, solving local ODEs of xi_immob
	INFO("--Solving local problem for xi_local and concentrations:");



	// loop over all the nodes
    for (node_idx = _concentrations[0]->getDiscreteData()->getRangeBegin();
	     node_idx < _concentrations[0]->getDiscreteData()->getRangeEnd();
		 node_idx++ )
	{

		// skip the boundary nodes
		if ( ! this->_solution->isBCNode(node_idx) )
		{

			// on each node, get the right start value
			// get the right set of eta and xi
			for (i=0; i < _n_eta; i++)
				loc_eta[i] = this->_eta[i]->getValue(node_idx);
			// fill in eta_immob
			for (i=0; i < _n_eta_bar; i++)
				loc_etabar[i] = this->_eta_bar[i]->getValue(node_idx);
			for (i=0; i < _n_xi_global; i++)
				loc_xi_global[i] = this->_xi_global[i]->getValue(node_idx);
			for (i=0; i < _n_xi_local; i++)
				loc_xi_local[i] = this->_xi_local[i]->getValue(node_idx);

			for (i=0; i < _n_Comp; i++)
				loc_conc[i] = this->_concentrations[i]->getValue(node_idx);

//			// take natural log of mobile and non mineral concentrations
//			loc_conc_non_Min = loc_conc.head(_I_mob + _I_NMin_bar);
//			loc_conc_Min 	 = loc_conc.tail(_I_min);


			// xi global constrains
			loc_XiSorpTilde  = loc_xi_global.head(_n_xi_Sorp_tilde);
			loc_XiMinTilde   = loc_xi_global.segment(_n_xi_Sorp_tilde, _n_xi_Min_tilde);
			loc_XiKin        = loc_xi_global.segment(_n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp + _n_xi_Min, _n_xi_Kin);

			loc_XiBarKin     = loc_xi_local.tail(_n_xi_Kin_bar);
			loc_XiBarKin_old = loc_xi_local_old_dt.tail(_n_xi_Kin_bar);

			// total mass constrain ie xi global and eta and etabar
			vec_tot_mass_constrain.segment(0,_n_eta) 																						    =  loc_eta;
			vec_tot_mass_constrain.segment(_n_eta,_n_xi_Sorp_tilde) 																		    =  loc_XiSorpTilde;
			vec_tot_mass_constrain.segment(_n_eta + _n_xi_Sorp_tilde,_n_xi_Min_tilde) 														    =  loc_XiMinTilde;
			vec_tot_mass_constrain.segment(_n_eta + _n_xi_Sorp_tilde + _n_xi_Min_tilde, _n_xi_Kin)											    =  loc_XiKin;
			vec_tot_mass_constrain.segment(_n_eta + _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Kin,_n_eta_bar) 								    =  loc_etabar;
			vec_tot_mass_constrain.segment(_n_eta + _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Kin + _n_eta_bar, _n_xi_Kin_bar)				    =  loc_XiBarKin_old;

			ogsChem::LocalVector ln_conc_Mob, ln_conc_NonMin_bar, conc_Mob, conc_NonMin_bar, conc_nonMin, ln_conc_nonMin, vec_conc_updated, conc_Min_bar;
			conc_Mob 		    = ogsChem::LocalVector::Zero(_I_mob);
			ln_conc_Mob 		= ogsChem::LocalVector::Zero(_I_mob);
			conc_NonMin_bar 	= ogsChem::LocalVector::Zero(_I_NMin_bar);
			ln_conc_NonMin_bar 	= ogsChem::LocalVector::Zero(_I_NMin_bar);
			conc_nonMin 		= ogsChem::LocalVector::Zero(_I_mob + _I_NMin_bar);
			ln_conc_nonMin 		= ogsChem::LocalVector::Zero(_I_mob + _I_NMin_bar);
			vec_conc_updated	= ogsChem::LocalVector::Zero(_n_Comp);
			conc_Min_bar	    = ogsChem::LocalVector::Zero(_I_min);

		    //conc_Min_bar        = vec_unknowns.segment(_I_mob + _I_NMin_bar, _I_min );
		    //Xi_Kin_bar   	    = vec_unknowns.segment(_I_mob + _I_NMin_bar + _I_min, _n_xi_Kin_bar );

			conc_nonMin = loc_conc.head(_I_mob + _I_NMin_bar );
			ln_conc_nonMin = loc_conc.head(_I_mob + _I_NMin_bar );
		    // convert the ln mobile conc to mobile conc
			pSolve->cal_ln_conc_vec((_I_mob + _I_NMin_bar), conc_nonMin, ln_conc_nonMin);
			loc_conc.head(_I_mob + _I_NMin_bar )  = ln_conc_nonMin;


			vec_unknowns.head (_n_Comp)		  = loc_conc;
			vec_unknowns.tail (_n_xi_Kin_bar) = loc_XiBarKin;

#ifdef _DEBUG
	// debugging--------------------------
//			std::cout << "======================================== \n";
//	std::cout << "vec_unknowns Vector BEFORE solving local problem: \n";
//	std::cout << node_idx << std::endl;
//	std::cout << vec_unknowns << std::endl;
//	std::cout << "======================================== \n";
	// end of debugging-------------------
#endif




			// solve the local problem
			pSolve->solve_LocalProblem_Newton_LineSearch( vec_unknowns, vec_tot_mass_constrain, node_idx , dt, iter_tol, rel_tol, max_iter);


			// re assembling xi local
			vec_conc   	      =  vec_unknowns.head (_n_Comp);
		    vec_XiBarKin      =  vec_unknowns.tail (_n_xi_Kin_bar);

		    // convert the ln mobile conc to mobile conc
		    ln_conc_Mob 	  =  vec_conc.head(_I_mob);
		    pSolve->cal_exp_conc_vec(_I_mob,ln_conc_Mob, conc_Mob);

		    loc_xi_mobile     = _mat_c_mob_2_xi_mob * conc_Mob;
		    vec_xi_mob        =  loc_xi_mobile.head(_n_xi_Mob);

		     // convert the log nonmineral conc to nonmineral conc
		    ln_conc_NonMin_bar  = vec_conc.segment(_I_mob, _I_NMin_bar);
		    pSolve->cal_exp_conc_vec(_I_NMin_bar, ln_conc_NonMin_bar, conc_NonMin_bar);
		    vec_conc.segment(_I_mob, _I_NMin_bar)  = conc_NonMin_bar;

		    loc_xi_immobile   = _mat_c_immob_2_xi_immob * vec_conc.tail(_I_NMin_bar + _n_xi_Min);
		    vec_XiSorpBar     = loc_xi_immobile.head(_n_xi_Sorp_bar);

		    vec_XiMinBar	  = vec_conc.tail(_n_xi_Min_bar);

		    loc_xi_local_new.head   (_n_xi_Mob) 								 = vec_xi_mob;
		    loc_xi_local_new.segment(_n_xi_Mob, _n_xi_Sorp_bar)				     = vec_XiSorpBar;
		    loc_xi_local_new.segment(_n_xi_Mob + _n_xi_Sorp_bar, _n_xi_Min_bar)  = vec_XiMinBar;
		    loc_xi_local_new.tail   (_n_xi_Kin_bar)							     = vec_XiBarKin;

		    vec_conc_updated.head(_I_mob)				   =  conc_Mob;
		    vec_conc_updated.segment(_I_mob, _I_NMin_bar)  =  conc_NonMin_bar;
		    conc_Min_bar								   =  vec_conc.tail(_I_min);
		    vec_conc_updated.tail(_I_min) 			       =  conc_Min_bar;

#ifdef _DEBUG
	// debugging--------------------------
	std::cout << "loc_xi_local_new Vector: \n";
	std::cout << loc_xi_local_new << std::endl;

	std::cout << "======================================== \n";

	std::cout << "vec_conc_updated Vector: \n";
	std::cout << vec_conc_updated << std::endl;

	// end of debugging-------------------
#endif


			// collect the xi_local_new
			for (i=0; i < _n_xi_local; i++)
				//_xi_local_new[i]->setValue(node_idx, loc_xi_local_new[i]);
				_xi_local[i]->setValue(node_idx, loc_xi_local_new[i]);

			// collect new concentration values
			for (i=0; i < _n_Comp; i++)
				_concentrations[i]->setValue(node_idx, vec_conc_updated[i]);
		} // end of if
        else
        {
            // if it is boundary nodes, then xi_local_new is equal to xi_local.
            for (i=0; i < _n_xi_local; i++)
				//_xi_local_new[i]->setValue(node_idx, _xi_local[i]->getValue( node_idx ) );
            	_xi_local[i]->setValue(node_idx, _xi_local[i]->getValue( node_idx ) );

            // if it is boundary nodes, then conc_glob_new is equal to conc_glob.
            for (i=0; i < _n_Comp; i++)
            	_concentrations[i]->setValue(node_idx, _concentrations[i]->getValue( node_idx ) );
        }

	}  // end of for node_idx

    delete pSolve;
}


template <class T1, class T2>
void FunctionReductConc<T1, T2>::GlobalResidualAssembler(const NumLib::TimeStep & delta_t, const SolutionLib::SolutionVector & u_cur_xiglob, SolutionLib::SolutionVector & residual_global)
{
	size_t i, node_idx, indx_tmp, nnodes;
    _n_xi_Sorp_bar_li  = _ReductionGIA->get_n_xi_Sorp_bar_li();
    _n_xi_Sorp_bar_ld  = _ReductionGIA->get_n_xi_Sorp_bar_ld();
    MyDiscreteSystem* dis = 0;
    nnodes = dis->getMesh()->getNumberOfNodes();
    // current xi global
	MathLib::LocalVector loc_cur_xi_global, loc_cur_xi_Sorp_tilde, loc_cur_xi_Min_tilde,
						 loc_cur_xi_Sorp, loc_cur_xi_Min, loc_cur_xi_Kin;
						 //global_cur_xi_Sorp_tilde, global_cur_xi_Min_tilde,
						 //global_cur_xi_Sorp, global_cur_xi_Min, global_cur_xi_Kin;
	// previous xi global
	MathLib::LocalVector loc_pre_xi_global, loc_pre_xi_Sorp_tilde, loc_pre_xi_Min_tilde, loc_pre_xi_Sorp, loc_pre_xi_Min, loc_pre_xi_Kin;
	                     //global_pre_xi_Sorp_tilde, global_pre_xi_Min_tilde, global_pre_xi_Sorp, global_pre_xi_Min, global_pre_xi_Kin;

	// current xi local
	MathLib::LocalVector loc_cur_xi_local, loc_xi_local, loc_cur_xi_Mob, loc_cur_xi_Sorp_bar, loc_cur_xi_Min_bar, loc_cur_xi_Kin_bar, loc_cur_xi_Sorp_bar_li, loc_cur_xi_Sorp_bar_ld;

	// residual vectors
	MathLib::LocalVector res43, res44, res45, res46, local_residual;//, residual_global_res43, residual_global_res44, residual_global_res45, residual_global_res46, global_vec_LHS_sorp, global_vec_RHS_sorp, global_vec_LHS_min, global_vec_RHS_min;

	// eta vectors
	MathLib::LocalVector loc_cur_eta, loc_cur_eta_bar;

	// rate vector
	MathLib::LocalVector vec_Rate, vec_Rate_45, vec_Rate_46; //, global_vec_Rate_45, global_vec_Rate_46;

	MathLib::LocalMatrix _mat_A1sorp, _mat_A2sorpli, _mat_A2sorpld, _mat_A1min, _mat_Ald;

	// initialize the local vector
	//current xi global
	loc_cur_xi_global       		= LocalVector::Zero( _n_xi_global );
	loc_cur_xi_Sorp_tilde			= LocalVector::Zero( _n_xi_Sorp_tilde );
	loc_cur_xi_Min_tilde			= LocalVector::Zero( _n_xi_Min_tilde );
	loc_cur_xi_Sorp					= LocalVector::Zero( _n_xi_Sorp);
	loc_cur_xi_Min					= LocalVector::Zero( _n_xi_Min);
	loc_cur_xi_Kin					= LocalVector::Zero( _n_xi_Kin);
	//previous xi global
	loc_pre_xi_global       		= LocalVector::Zero( _n_xi_global );
	loc_pre_xi_Sorp_tilde			= LocalVector::Zero( _n_xi_Sorp_tilde );
	loc_pre_xi_Min_tilde			= LocalVector::Zero( _n_xi_Min_tilde );
	loc_pre_xi_Sorp					= LocalVector::Zero( _n_xi_Sorp);
	loc_pre_xi_Min					= LocalVector::Zero( _n_xi_Min);
	loc_pre_xi_Kin					= LocalVector::Zero( _n_xi_Kin);
	//current xi local
	loc_cur_xi_local       			= LocalVector::Zero( _n_xi_local );
	loc_cur_xi_Sorp_bar        	    = LocalVector::Zero( _n_xi_Sorp_bar );
	loc_cur_xi_Min_bar       		= LocalVector::Zero( _n_xi_Min_bar );
	loc_cur_xi_Sorp_bar_li			= LocalVector::Zero(_n_xi_Sorp_bar_li );
	loc_cur_xi_Sorp_bar_ld			= LocalVector::Zero(_n_xi_Sorp_bar_ld );
	//residual vectors
	res43						= LocalVector::Zero( _n_xi_Sorp_tilde);
	res44						= LocalVector::Zero( _n_xi_Min_tilde);
	res45						= LocalVector::Zero( _n_xi_Sorp);
	res46						= LocalVector::Zero( _n_xi_Min);
//	residual_global_res43		= LocalVector::Zero( _n_xi_Sorp_tilde * nnodes);
//	residual_global_res44		= LocalVector::Zero( _n_xi_Min_tilde * nnodes);
//	residual_global_res45		= LocalVector::Zero( _n_xi_Min * nnodes);
//	residual_global_res46		= LocalVector::Zero( _n_xi_Sorp * nnodes);
//	global_vec_LHS_sorp			= LocalVector::Zero( _n_xi_Sorp * nnodes);
//	global_vec_RHS_sorp			= LocalVector::Zero( _n_xi_Sorp * nnodes);
//	global_vec_LHS_min			= LocalVector::Zero( _n_xi_Min * nnodes);
//	global_vec_RHS_min			= LocalVector::Zero( _n_xi_Min * nnodes);
	// current eta mobie and immobile
	loc_cur_eta        		= LocalVector::Zero( _n_eta );
	loc_cur_eta_bar     	= LocalVector::Zero( _n_eta_bar );


	//initialize it before?
    _mat_A1sorp		   = _ReductionGIA->get_matrix_A1sorp();
    _mat_A2sorpli 	   = _ReductionGIA->get_matrix_A2sorpli();
    _mat_A1min		   = _ReductionGIA->get_matrix_A1min();
    _mat_Ald 		   = _ReductionGIA->get_matrix_Ald();
    _mat_A2sorpld  	   = _ReductionGIA->get_matrix_A2sorpld();

    _mat_Asorp = _mat_A1sorp - _mat_A2sorpli;
    _mat_Amin = _mat_A1min - _mat_Ald * _mat_A2sorpld;


	// loop over all the nodes
    for (node_idx = _concentrations[0]->getDiscreteData()->getRangeBegin();
	     node_idx < _concentrations[0]->getDiscreteData()->getRangeEnd();
		 node_idx++ )
	{
		// skip the boundary nodes
		if ( ! this->_solution->isBCNode(node_idx) )
		{

			// on each node, get the right start value
			// get the right set of xis

			for (i=0; i < _n_xi_global; i++)
				loc_cur_xi_global[i] = u_cur_xiglob[node_idx * _n_xi_global + i];
			for (i=0; i < _n_xi_global; i++)
				loc_pre_xi_global[i] = this->_xi_global[i]->getValue(node_idx);
			for (i=0; i < _n_xi_local; i++)
				loc_cur_xi_local[i] = this->_xi_local[i]->getValue(node_idx);
            for (i=0; i < _n_eta; i++)
                loc_cur_eta[i] = this->_eta[i]->getValue(node_idx);
            // fill in eta_immob
            for (i=0; i < _n_eta_bar; i++)
                loc_cur_eta_bar[i] = this->_eta_bar[i]->getValue(node_idx);


		    // current xi global
			loc_cur_xi_Sorp_tilde   = loc_cur_xi_global.head(_n_xi_Sorp_tilde);
			loc_cur_xi_Min_tilde  	= loc_cur_xi_global.segment(_n_xi_Sorp_tilde, _n_xi_Min_tilde);
		    loc_cur_xi_Sorp  	 	= loc_cur_xi_global.segment( _n_xi_Sorp_tilde + _n_xi_Min_tilde, _n_xi_Sorp);
		    loc_cur_xi_Min   		= loc_cur_xi_global.segment( _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp, _n_xi_Min);
		    loc_cur_xi_Kin  		= loc_cur_xi_global.tail(_n_xi_Kin);


//			for (i=0; i < _n_xi_Sorp_tilde; i++)
//				_global_cur_xi_Sorp_tilde[i]->setValue(node_idx, loc_cur_xi_Sorp_tilde[i]);
//			for (i=0; i < _n_xi_Min_tilde; i++)
//				_global_cur_xi_Min_tilde[i]->setValue(node_idx, loc_cur_xi_Min_tilde[i]);
//			for (i=0; i < _n_xi_Sorp; i++)
//				_global_cur_xi_Sorp[i]->setValue(node_idx, loc_cur_xi_Sorp[i]);
//			for (i=0; i < _n_xi_Min; i++)
//				_global_cur_xi_Min[i]->setValue(node_idx, loc_cur_xi_Min[i]);
//			for (i=0; i < _n_xi_Kin; i++)
//				_global_cur_xi_Kin[i]->setValue(node_idx, loc_cur_xi_Kin[i]);


			// previous xi global
			loc_pre_xi_Sorp_tilde 	= loc_pre_xi_global.head(_n_xi_Sorp_tilde);
			loc_pre_xi_Min_tilde  	= loc_pre_xi_global.segment(_n_xi_Sorp_tilde, _n_xi_Min_tilde);
		    loc_pre_xi_Sorp   		= loc_pre_xi_global.segment( _n_xi_Sorp_tilde + _n_xi_Min_tilde, _n_xi_Sorp);
		    loc_pre_xi_Min   		= loc_pre_xi_global.segment( _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp, _n_xi_Min);
		    loc_pre_xi_Kin  		= loc_pre_xi_global.tail(_n_xi_Kin);

//			for (i=0; i < _n_xi_Sorp_tilde; i++)
//				_global_pre_xi_Sorp_tilde[i]->setValue(node_idx, loc_pre_xi_Sorp_tilde[i]);
//			for (i=0; i < _n_xi_Min_tilde; i++)
//				_global_pre_xi_Min_tilde[i]->setValue(node_idx, loc_pre_xi_Min_tilde[i]);
//			for (i=0; i < _n_xi_Sorp; i++)
//				_global_pre_xi_Sorp[i]->setValue(node_idx, loc_pre_xi_Sorp[i]);
//			for (i=0; i < _n_xi_Min; i++)
//				_global_pre_xi_Min[i]->setValue(node_idx, loc_pre_xi_Min[i]);
//			for (i=0; i < _n_xi_Kin; i++)
//				_global_pre_xi_Kin[i]->setValue(node_idx, loc_pre_xi_Kin[i]);


			// current xi local
			loc_cur_xi_Mob   		= loc_cur_xi_local.head(_n_xi_Mob);
			loc_cur_xi_Sorp_bar 	= loc_cur_xi_local.segment(_n_xi_Mob, _n_xi_Sorp_bar);
			loc_cur_xi_Min_bar 		= loc_cur_xi_local.segment( _n_xi_Mob + _n_xi_Sorp_bar, _n_xi_Min_bar);
			loc_cur_xi_Kin_bar 		= loc_cur_xi_local.tail(_n_xi_Kin_bar);

		    loc_cur_xi_Sorp_bar_li 	= loc_cur_xi_Sorp_bar.topRows(_n_xi_Sorp_bar_li);
		    loc_cur_xi_Sorp_bar_ld 	= loc_cur_xi_Sorp_bar.bottomRows(_n_xi_Sorp_bar_ld);





		    // calculate node based AE (Eq. 3.43 and 3.44)
		    res43 = loc_cur_xi_Sorp_tilde - loc_cur_xi_Sorp + loc_cur_xi_Sorp_bar_li;
		    for(std::size_t i = 0; i < _n_xi_Sorp_tilde; i++)
		    residual_global[_n_xi_global * node_idx + i] = res43[i];

		    res44 = loc_cur_xi_Min_tilde - loc_cur_xi_Min + loc_cur_xi_Min_bar + loc_cur_xi_Sorp_bar_ld;
		    for(std::size_t i = 0; i < _n_xi_Min_tilde; i++)
		    residual_global[_n_xi_global * node_idx + i] = res44[i];

//			// collect the xi_local_new  ??
//			for (i=0; i < _n_xi_Sorp_tilde; i++)
//				residual_global_res43[i]->setValue(node_idx, res43[i]);
//			// residual_global_res43.push_back( res43 );
//
//			for (i=0; i < _n_xi_Min_tilde; i++)
//				residual_global_res44[i]->setValue(node_idx, res44[i]);
//			// residual_global_res44.push_back( res44 );


			// Note: since rate is already used in local problem, is it the same? probably not.
			// calculate the nodal kinetic reaction rates
			this->_ReductionGIA->Calc_Kin_Rate(loc_cur_xi_Mob,
											loc_pre_xi_Sorp,
											loc_cur_xi_Sorp_tilde,
											loc_cur_xi_Sorp_bar,
											loc_cur_xi_Min,
											loc_cur_xi_Min_tilde,
											loc_cur_xi_Min_bar,
											loc_cur_xi_Kin,
											loc_cur_xi_Kin_bar,
											loc_cur_eta,
											loc_cur_eta_bar,
											vec_Rate);



			_vec_Rate_rows  = vec_Rate.rows();
			for (i=0; i < vec_Rate.rows(); i++)
				_global_vec_Rate[i]->setValue(node_idx, vec_Rate[i]);

//			for (i=0; i < vec_Rate_46.cols(); i++)
//				_global_vec_Rate_46[i]->setValue(node_idx, vec_Rate_46[i]);

		}

	}

    // solve for Eq. 45
    std::size_t switchOn = 1;  // 1 for xi sorp tilde
    this->assembly(delta_t, _n_xi_Sorp, u_cur_xiglob, residual_global, switchOn);

    switchOn = 0;  // 1 for xi sorp tilde
    this->assembly(delta_t, _n_xi_Min, u_cur_xiglob, residual_global, switchOn);

//    for(std::vector<int>::iterator iter; iter _residual_global_res45.begin(); iter != _residual_global_res45.end(); ++iter)
//    _residual_global_res45[iter] = _global_vec_LHS_sorp[iter] - _global_vec_RHS_sorp[iter] - _global_vec_Rate_45[iter];
//
//    // solve for Eq. 46
//    this->assembly(delta_t, _n_xi_Min, global_cur_xi_Min_tilde, global_pre_xi_Min_tilde, global_cur_xi_Min, global_pre_xi_Min, _global_vec_LHS_min, _global_vec_RHS_min);
//
//    residual_global_res46 = _global_vec_LHS_min - _global_vec_RHS_min - _global_vec_Rate_46;
//
//
//
//    // construct global residual vector
//	// loop over all the nodes
//    for (node_idx = _concentrations[0]->getDiscreteData()->getRangeBegin();
//	     node_idx < _concentrations[0]->getDiscreteData()->getRangeEnd();
//		 node_idx++ )
//	{
//
//	for (i=0; i < _n_xi_Sorp_tilde; i++)
//		res43[i] = this->_residual_global_res43[i]->getValue(node_idx);
//	for (i=0; i < _n_xi_Min_tilde; i++)
//		res44[i] = this->_residual_global_res44[i]->getValue(node_idx);
//
//	for (i=0; i < _n_xi_Sorp; i++)
//		res45[i] = this->_residual_global_res45[i]->getValue(node_idx);
//	for (i=0; i < _n_xi_Min; i++)
//		res46[i] = this->_residual_global_res46[i]->getValue(node_idx);
//
//	local_residual.head(_n_xi_Sorp_tilde) 									= res43;
//	local_residual.segment(_n_xi_Sorp_tilde,_n_xi_Min_tilde) 				= res44;
//	local_residual.segment(_n_xi_Sorp_tilde + _n_xi_Min_tilde, _n_xi_Sorp) 	= res45;
//	local_residual.tail(_n_xi_Min) 											= res46;
//
//	for (i=0; i < _n_xi_global; i++)
//		_residual_global[i]->setValue(node_idx, local_residual[i]);
//
//
//	}
}


// element based assembly
template <class T1, class T2>
void FunctionReductConc<T1, T2>::assembly(const NumLib::TimeStep & delta_t, const std::size_t _n_xi, const SolutionLib::SolutionVector & u_cur_xiglob,
											SolutionLib::SolutionVector & residual_global, std::size_t switchOn)
		{
	//MyLinearSolver* linear_solver = this->_non_linear_solution->getLinearEquationSolver();
	MeshLib::IMesh* msh = _dis_sys->getMesh();
	const size_t n_ele = msh->getNumberOfElements();
	double _theta(1.0);

	MathLib::LocalVector loc_cur_xi_global, loc_cur_xi_Sorp_tilde, loc_cur_xi_Min_tilde, loc_cur_xi_Sorp, loc_cur_xi_Min, loc_cur_xi_Kin;

	// previous xi global
	MathLib::LocalVector loc_pre_xi_global, loc_pre_xi_Sorp_tilde, loc_pre_xi_Min_tilde, loc_pre_xi_Sorp, loc_pre_xi_Min, loc_pre_xi_Kin;

	MathLib::LocalVector loc_vec_Rate, vec_Rate_temp;

	// initialize the local vector
	//current xi global
	loc_cur_xi_global       		= LocalVector::Zero( _n_xi_global );
	loc_cur_xi_Sorp_tilde			= LocalVector::Zero( _n_xi_Sorp_tilde );
	loc_cur_xi_Min_tilde			= LocalVector::Zero( _n_xi_Min_tilde );
	loc_cur_xi_Sorp					= LocalVector::Zero( _n_xi_Sorp);
	loc_cur_xi_Min					= LocalVector::Zero( _n_xi_Min);
	loc_cur_xi_Kin					= LocalVector::Zero( _n_xi_Kin);
	//previous xi global
	loc_pre_xi_global       		= LocalVector::Zero( _n_xi_global );
	loc_pre_xi_Sorp_tilde			= LocalVector::Zero( _n_xi_Sorp_tilde );
	loc_pre_xi_Min_tilde			= LocalVector::Zero( _n_xi_Min_tilde );
	loc_pre_xi_Sorp					= LocalVector::Zero( _n_xi_Sorp);
	loc_pre_xi_Min					= LocalVector::Zero( _n_xi_Min);
	loc_pre_xi_Kin					= LocalVector::Zero( _n_xi_Kin);


	for (size_t i=0; i<n_ele; i++)
	{
		MeshLib::IElement *e = msh->getElement(i);

	std::vector<size_t> ele_node_ids;
	//MathLib::Vector<size_t> ele_node_ids;
	e->getNodeIDList(e->getMaximumOrder(), ele_node_ids);

	//TODO get the water content and multiply it in LHS and RHS

	double dt = delta_t.getTimeStepSize();

	MathLib::LocalMatrix localM = MathLib::LocalMatrix::Zero(ele_node_ids.size(), ele_node_ids.size());
	MathLib::LocalMatrix localK = MathLib::LocalMatrix::Zero(ele_node_ids.size(), ele_node_ids.size());
	MathLib::LocalVector F = MathLib::LocalVector::Zero(ele_node_ids.size());

	//assembleODE(time, e, local_u_n1, local_u_n, M, K, F);
	//_fe = _feObjects->getFeObject(e);

	const size_t n_dim = e->getDimension();
	size_t mat_id = e->getGroupID();
	_pm = Ogs6FemData::getInstance()->list_pm[mat_id];

	localDispersion.setZero(localK.rows(), localK.cols());
	localAdvection.setZero (localK.rows(), localK.cols());

	double cmp_mol_diffusion = .0;
	// _cmp->molecular_diffusion->eval(0, cmp_mol_diffusion);

	_q = _fe->getIntegrationMethod();
	double gp_x[3], real_x[3];
	poro.setZero();
	d_poro.setZero();
	double disp_l = 0.0;
	double disp_t = 0.0;

	for (size_t j=0; j < _q->getNumberOfSamplingPoints(); j++)
	{
		_q->getSamplingPoint(j, gp_x);
		_fe->computeBasisFunctions(gp_x);
		_fe->getRealCoordinates(real_x);
		NumLib::TXPosition gp_pos(NumLib::TXPosition::IntegrationPoint, e->getID(), j, real_x);

		_pm->porosity->eval(gp_pos, poro);
		_pm->dispersivity_long->eval(gp_pos, disp_l);
		_pm->dispersivity_trans->eval(gp_pos, disp_t);

		d_poro(0,0) = cmp_mol_diffusion * poro(0,0);
		d_poro(1,1) = cmp_mol_diffusion * poro(0,0);
		d_poro(2,2) = cmp_mol_diffusion * poro(0,0);
		_vel->eval(gp_pos, v);
		v2 = v.topRows(n_dim).transpose();

		// calculating dispersion tensor according to Benchmark book p219, Eq. 10.15
		// D_{ij} = \alpha_T |v| \delta_{ij} + (\alpha_L - \alpha_T) \frac{v_i v_j}{|v|} + D^{d}_{ii}
		dispersion_diffusion.setIdentity(n_dim, n_dim);
		dispersion_diffusion *= disp_l * v.norm();
		dispersion_diffusion += (disp_l - disp_t) * ( v2.transpose() * v2 ) / v.norm();
		dispersion_diffusion += d_poro.topLeftCorner(n_dim, n_dim);
		// --------debugging--------------
		// std::cout << "dispersion_diffusion Matrix" << std::endl;
		// std::cout << dispersion_diffusion << std::endl;
		// --------end of debugging-------

		_fe->integrateWxN(j, poro, localM);
		_fe->integrateDWxDN(j, dispersion_diffusion, localDispersion);
		_fe->integrateWxDN(j, v2, localAdvection);
	}

	localK = localDispersion + localAdvection;

	// mass lumping----------------------------
	for (size_t idx_ml=0; idx_ml < localM.rows(); idx_ml++ )
	{
		double mass_lump_val;
		mass_lump_val = localM.row(idx_ml).sum();
		localM.row(idx_ml).setZero();
		localM(idx_ml, idx_ml) = mass_lump_val;
	}

	//std::cout << "M="; M.write(std::cout); std::cout << std::endl;
	//std::cout << "K="; K.write(std::cout); std::cout << std::endl;

	//MathLib::LocalVector localLHS = LocalVector::Zero(ele_node_ids.size());
	//MathLib::LocalVector localRHS = LocalVector::Zero(ele_node_ids.size());
	MathLib::LocalVector local_u1 = LocalVector::Zero(ele_node_ids.size());
	MathLib::LocalVector local_u0 = LocalVector::Zero(ele_node_ids.size());
	MathLib::LocalVector local_u1_tilde = LocalVector::Zero(ele_node_ids.size());
	MathLib::LocalVector local_u0_tilde = LocalVector::Zero(ele_node_ids.size());

	//MathLib::LocalVector node_indx = LocalVector::Zero(_n_xi);
	std::size_t xi_count, node_indx(0);
	double localLHS, localRHS;

	for (xi_count = 0; xi_count < _n_xi; xi_count++)
	{

//		for(size_t idx = 0; idx != ele_node_ids.size(); idx++)
//		{
//
//			node_indx			= ele_node_ids[idx] * _n_xi + xi_count;
//			local_u1(idx) 		= _global_u1_tilde[node_indx];
//			local_u0(idx) 		= _global_u0[node_indx];
//			local_u1_tilde(idx) = _global_u1_tilde[node_indx];
//			local_u0_tilde(idx) = _global_u0_tilde[node_indx];
//		}


		for(size_t idx = 0; idx != ele_node_ids.size(); idx++)
		{
			node_indx = ele_node_ids[idx] * _n_xi + xi_count;

			for (i=0; i < _n_xi_global; i++)
				loc_cur_xi_global[i] = u_cur_xiglob[ele_node_ids[idx] * _n_xi_global + xi_count];
			for (i=0; i < _n_xi_global; i++)
				loc_pre_xi_global[i] = this->_xi_global[i]->getValue(ele_node_ids[idx]);

			loc_cur_xi_Sorp_tilde   = loc_cur_xi_global.head(_n_xi_Sorp_tilde);
			loc_cur_xi_Min_tilde  	= loc_cur_xi_global.segment(_n_xi_Sorp_tilde, _n_xi_Min_tilde);
		    loc_cur_xi_Sorp  	 	= loc_cur_xi_global.segment( _n_xi_Sorp_tilde + _n_xi_Min_tilde, _n_xi_Sorp);
		    loc_cur_xi_Min   		= loc_cur_xi_global.segment( _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp, _n_xi_Min);

			loc_pre_xi_Sorp_tilde 	= loc_pre_xi_global.head(_n_xi_Sorp_tilde);
			loc_pre_xi_Min_tilde  	= loc_pre_xi_global.segment(_n_xi_Sorp_tilde, _n_xi_Min_tilde);
		    loc_pre_xi_Sorp   		= loc_pre_xi_global.segment( _n_xi_Sorp_tilde + _n_xi_Min_tilde, _n_xi_Sorp);
		    loc_pre_xi_Min   		= loc_pre_xi_global.segment( _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp, _n_xi_Min);

			for (i=0; i < _vec_Rate_rows; i++)
				loc_vec_Rate[i] = this->_global_vec_Rate[i]->getValue(ele_node_ids[idx]);

			if(switchOn)  //if switchOn calculate for xi sorp tilde and for xi min tilde otherwise
			{
							local_u1(idx) 		= loc_cur_xi_Sorp[node_indx];
							local_u0(idx) 		= loc_pre_xi_Sorp[node_indx];
							local_u1_tilde(idx) = loc_cur_xi_Sorp_tilde[node_indx];
							local_u0_tilde(idx) = loc_pre_xi_Sorp_tilde[node_indx];
							vec_Rate_temp = _mat_Asorp * loc_vec_Rate;

			}
			else
			{
							local_u1(idx) 		= loc_cur_xi_Min[node_indx];
							local_u0(idx) 		= loc_pre_xi_Min[node_indx];
							local_u1_tilde(idx) = loc_cur_xi_Min_tilde[node_indx];
							local_u0_tilde(idx) = loc_pre_xi_Min_tilde[node_indx];
							vec_Rate_temp = _mat_Amin * loc_vec_Rate;
			}

		// LHS = dt*(1/dt M + theta K)
		//localLHS = (localM * global_u1_tilde(node_indx)) + (delta_t * _theta * localK * global_u1(node_indx));
		 localLHS = (localM.row(idx).dot(local_u1_tilde)) + (dt * _theta * (localK.row(idx).dot(local_u1)));
		//_global_vec_LHS(node_indx) = _global_vec_LHS(node_indx) + localLHS;
		 // RHS = (1/dt M - (1-theta) K) u0 + F
		 localRHS = (localM.row(idx).dot(local_u0_tilde)) - (dt * (1.-_theta) * localK.row(idx).dot(local_u0)); //local_u_n;
		 //_global_vec_RHS(node_indx) = _global_vec_RHS(node_indx) + localRHS;

		residual_global[_n_xi_global * ele_node_ids[idx] + xi_count] = localLHS - localRHS - vec_Rate_temp[xi_count];

		}
	}
	}
		}
//
//	 Map<RowVectorXi> local_u1(global_u1,node_indx);
//	 Map<RowVectorXi> local_u0(global_u0,node_indx);
//	 Map<RowVectorXi> local_u1_tilde(global_u1_tilde,node_indx);
//	 Map<RowVectorXi> local_u0_tilde(global_u0_tilde,node_indx);
//
//	// LHS = dt*(1/dt M + theta K)
//	//localLHS = (localM * global_u1_tilde(node_indx)) + (delta_t * _theta * localK * global_u1(node_indx));
//	 localLHS = (localM * local_u1_tilde) + (delta_t * _theta * localK * local_u1);
//	_global_vec_LHS(node_indx) = _global_vec_LHS(node_indx) + localLHS;
//	// RHS = (1/dt M - (1-theta) K) u0 + F
//	localRHS = (localM * local_u0_tilde) - (delta_t * (1.-_theta) * localK * local_u0); //local_u_n;
//	_global_vec_RHS(node_indx) = _global_vec_RHS(node_indx) + localRHS;
//	}
//
//
//	}
//		}






template <class T1, class T2>
void FunctionReductConc<T1, T2>::GlobalJacobianAssembler(const NumLib::TimeStep & delta_t, const SolutionLib::SolutionVector & u_cur_xiglob, SolutionLib::SolutionVector & Jacobian_global)  // TODO Jacobian_global will be changed to matrix
{

	using namespace std::placeholders;
	size_t i, node_idx, indx_tmp, nnodes;
    _n_xi_Sorp_bar_li  = _ReductionGIA->get_n_xi_Sorp_bar_li();
    _n_xi_Sorp_bar_ld  = _ReductionGIA->get_n_xi_Sorp_bar_ld();
    MyDiscreteSystem* dis = 0;
    nnodes = dis->getMesh()->getNumberOfNodes();
    const double theta_water_content(0.32);
    const double delta_xi = 1E-14;

    MathLib::LocalMatrix mat_A1kin   = _ReductionGIA->get_matrix_A1kin();
    ogsChem::LocalMatrix mat_Ald   = _ReductionGIA->get_matrix_Ald();

    // current xi global
	MathLib::LocalVector loc_cur_xi_global, loc_cur_xi_Sorp_tilde, loc_cur_xi_Min_tilde,
						 loc_cur_xi_Sorp, loc_cur_xi_Min, loc_cur_xi_Kin;
	// previous xi global
	MathLib::LocalVector loc_pre_xi_global, loc_pre_xi_Sorp_tilde, loc_pre_xi_Min_tilde, loc_pre_xi_Sorp, loc_pre_xi_Min, loc_pre_xi_Kin;

	// current xi local
	MathLib::LocalVector loc_cur_xi_local, loc_xi_local, loc_cur_xi_Mob, loc_cur_xi_Sorp_bar, loc_cur_xi_Min_bar, loc_cur_xi_Kin_bar, loc_cur_xi_Sorp_bar_li, loc_cur_xi_Sorp_bar_ld;

	// eta vectors
	MathLib::LocalVector loc_cur_eta, loc_cur_eta_bar;

	// local concentration vector
	MathLib::LocalVector vec_conc;

	MathLib::LocalMatrix mat_p1Fder,mat_p1Ftrans, mat_p1F, mat_p2F, mat_Global_Jacobian;

	// initialize the local vector
	//current xi global
	loc_cur_xi_global       		= LocalVector::Zero( _n_xi_global );
	loc_cur_xi_Sorp_tilde			= LocalVector::Zero( _n_xi_Sorp_tilde );
	loc_cur_xi_Min_tilde			= LocalVector::Zero( _n_xi_Min_tilde );
	loc_cur_xi_Sorp					= LocalVector::Zero( _n_xi_Sorp);
	loc_cur_xi_Min					= LocalVector::Zero( _n_xi_Min);
	loc_cur_xi_Kin					= LocalVector::Zero( _n_xi_Kin);
	//current xi local
	loc_cur_xi_local       			= LocalVector::Zero( _n_xi_local );
	loc_cur_xi_Sorp_bar        	    = LocalVector::Zero( _n_xi_Sorp_bar );
	loc_cur_xi_Min_bar       		= LocalVector::Zero( _n_xi_Min_bar );
	loc_cur_xi_Sorp_bar_li			= LocalVector::Zero(_n_xi_Sorp_bar_li );
	loc_cur_xi_Sorp_bar_ld			= LocalVector::Zero(_n_xi_Sorp_bar_ld );
	// current eta mobie and immobile
	loc_cur_eta        		= LocalVector::Zero( _n_eta );
	loc_cur_eta_bar     	= LocalVector::Zero( _n_eta_bar );
	vec_conc                = LocalVector::Zero(_n_Comp);

	MathLib::LocalMatrix mat_vprime  = LocalMatrix::Zero(_n_xi_local,_n_xi_global);
	MathLib::LocalMatrix Jacobian_local = LocalMatrix::Zero(_n_xi_global, _n_xi_global);
	MathLib::LocalVector vec_rate_old = LocalVector::Zero(_vec_Rate_rows);
	MathLib::LocalVector vec_Rate = LocalVector::Zero(_vec_Rate_rows);

    mat_p1Fder 	= LocalMatrix::Zero(_n_xi_global, _n_xi_global);
    mat_p1Ftrans = LocalMatrix::Zero(_n_xi_global, _n_xi_global);
    mat_p1F		= LocalMatrix::Zero(_n_xi_global, _n_xi_global);
    mat_p2F 	= LocalMatrix::Zero(_n_xi_global, _n_xi_global);
    double dt = delta_t.getTimeStepSize();

    //TODO fix the global jacobian matrix

	// loop over all the nodes
    for (node_idx = _concentrations[0]->getDiscreteData()->getRangeBegin();
	     node_idx < _concentrations[0]->getDiscreteData()->getRangeEnd();
		 node_idx++ )
	{
		// skip the boundary nodes
		if ( ! this->_solution->isBCNode(node_idx) )
		{

			// on each node, get the right start value
			// get the right set of xis

			for (i=0; i < _n_xi_global; i++)
				loc_cur_xi_global[i] = u_cur_xiglob[node_idx * _n_xi_global + i];
			for (i=0; i < _n_xi_local; i++)
				loc_cur_xi_local[i] = this->_xi_local[i]->getValue(node_idx);
            for (i=0; i < _n_eta; i++)
                loc_cur_eta[i] = this->_eta[i]->getValue(node_idx);
            // fill in eta_immob
            for (i=0; i < _n_eta_bar; i++)
                loc_cur_eta_bar[i] = this->_eta_bar[i]->getValue(node_idx);

			for (i=0; i < _vec_Rate_rows; i++)
				vec_rate_old[i] = this->_global_vec_Rate[i]->getValue(node_idx);

			for (i=0; i < _n_Comp; i++)
				vec_conc[i] = this->_concentrations[i]->getValue(node_idx);

		    // current xi global
			loc_cur_xi_Sorp_tilde   = loc_cur_xi_global.head(_n_xi_Sorp_tilde);
			loc_cur_xi_Min_tilde  	= loc_cur_xi_global.segment(_n_xi_Sorp_tilde, _n_xi_Min_tilde);
		    loc_cur_xi_Sorp  	 	= loc_cur_xi_global.segment( _n_xi_Sorp_tilde + _n_xi_Min_tilde, _n_xi_Sorp);
		    loc_cur_xi_Min   		= loc_cur_xi_global.segment( _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp, _n_xi_Min);
		    loc_cur_xi_Kin  		= loc_cur_xi_global.tail(_n_xi_Kin);

			// current xi local
			loc_cur_xi_Mob   		= loc_cur_xi_local.head(_n_xi_Mob);
			loc_cur_xi_Sorp_bar 	= loc_cur_xi_local.segment(_n_xi_Mob, _n_xi_Sorp_bar);
			loc_cur_xi_Min_bar 		= loc_cur_xi_local.segment( _n_xi_Mob + _n_xi_Sorp_bar, _n_xi_Min_bar);
			loc_cur_xi_Kin_bar 		= loc_cur_xi_local.tail(_n_xi_Kin_bar);

		    loc_cur_xi_Sorp_bar_li 	= loc_cur_xi_Sorp_bar.topRows(_n_xi_Sorp_bar_li);
		    loc_cur_xi_Sorp_bar_ld 	= loc_cur_xi_Sorp_bar.bottomRows(_n_xi_Sorp_bar_ld);


		    ///// calculate partial1F

			ogsChem::LocalVector der_sorpT_R;
			der_sorpT_R  = LocalVector::Zero(_n_xi_Sorp_tilde);
			NumDiff(_n_xi_Sorp_tilde, delta_xi, std::bind(this->_ReductionGIA->Calc_Kin_Rate ,loc_cur_xi_Mob, loc_cur_xi_Sorp, _3, loc_cur_xi_Sorp_bar, loc_cur_xi_Min, loc_cur_xi_Min_tilde,
					 loc_cur_xi_Min_bar, loc_cur_xi_Kin, loc_cur_xi_Kin_bar, loc_cur_eta, loc_cur_eta_bar, vec_Rate), vec_rate_old, loc_cur_xi_Sorp_tilde, der_sorpT_R);


			ogsChem::LocalVector der_minT_R;
			der_minT_R  = LocalVector::Zero(_n_xi_Min_tilde);
			NumDiff(_n_xi_Min_tilde, delta_xi, std::bind(this->_ReductionGIA->Calc_Kin_Rate ,loc_cur_xi_Mob, loc_cur_xi_Sorp, loc_cur_xi_Sorp_tilde, loc_cur_xi_Sorp_bar, loc_cur_xi_Min, _6,
					 loc_cur_xi_Min_bar, loc_cur_xi_Kin, loc_cur_xi_Kin_bar, loc_cur_eta, loc_cur_eta_bar, vec_Rate), vec_rate_old, loc_cur_xi_Min_tilde, der_minT_R);

			ogsChem::LocalVector der_kin_R;
			der_kin_R  = LocalVector::Zero(_n_xi_Kin);
			NumDiff(_n_xi_Kin, delta_xi, std::bind(this->_ReductionGIA->Calc_Kin_Rate ,loc_cur_xi_Mob, loc_cur_xi_Sorp, loc_cur_xi_Sorp_tilde, loc_cur_xi_Sorp_bar, loc_cur_xi_Min, loc_cur_xi_Min_tilde,
					 loc_cur_xi_Min_bar, _8, loc_cur_xi_Kin_bar, loc_cur_eta, loc_cur_eta_bar, vec_Rate), vec_rate_old, loc_cur_xi_Kin, der_kin_R);




			mat_p1Fder.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde, 0, _n_xi_Sorp, _n_xi_Sorp_tilde)						  =  _mat_Asorp * der_sorpT_R;
			mat_p1Fder.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp, 0, _n_xi_Min, _n_xi_Sorp_tilde) 			  =  _mat_Amin * der_sorpT_R;
			mat_p1Fder.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp + _n_xi_Min, 0, _n_xi_Kin, _n_xi_Sorp_tilde) =  mat_A1kin * der_sorpT_R;

			mat_p1Fder.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde, _n_xi_Sorp_tilde, _n_xi_Sorp, _n_xi_Min_tilde) 						=  _mat_Asorp * der_minT_R;
			mat_p1Fder.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp, _n_xi_Sorp_tilde, _n_xi_Min, _n_xi_Min_tilde) 			=  _mat_Amin * der_minT_R;
			mat_p1Fder.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp + _n_xi_Min, _n_xi_Sorp_tilde, _n_xi_Kin, _n_xi_Min_tilde) =  mat_A1kin * der_minT_R;

			mat_p1Fder.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde, _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp + _n_xi_Min, _n_xi_Sorp, _n_xi_Kin) 						 =  _mat_Asorp * der_kin_R;
			mat_p1Fder.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp, _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp + _n_xi_Min, _n_xi_Min, _n_xi_Kin) 			 =  _mat_Amin * der_kin_R;
			mat_p1Fder.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp + _n_xi_Min, _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp + _n_xi_Min, _n_xi_Kin, _n_xi_Kin) =  mat_A1kin * der_kin_R;

		    ///// Add the identity matrix to the  the mass and conductance matrix
			mat_p1Ftrans.block(0, 0, _n_xi_Sorp_tilde, _n_xi_Sorp_tilde) 													  = LocalMatrix::Identity(_n_xi_Sorp_tilde, _n_xi_Sorp_tilde);
			mat_p1Ftrans.block(_n_xi_Sorp_tilde, _n_xi_Sorp_tilde, _n_xi_Min_tilde, _n_xi_Min_tilde) 						  = LocalMatrix::Identity(_n_xi_Min_tilde, _n_xi_Min_tilde);
			mat_p1Ftrans.block(0, _n_xi_Sorp_tilde + _n_xi_Min_tilde, _n_xi_Sorp_tilde, _n_xi_Sorp)							  = -1.0 * LocalMatrix::Identity(_n_xi_Sorp_tilde, _n_xi_Sorp);
			mat_p1Ftrans.block(_n_xi_Sorp_tilde, _n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp, _n_xi_Min_tilde, _n_xi_Min) = -1.0 * LocalMatrix::Identity(_n_xi_Min_tilde, _n_xi_Min);

			//// calculate partial1f. i.e. partial drivative with respect to global variable
			mat_p1F = mat_p1Ftrans - dt * theta_water_content * mat_p1Fder;

			///// calculate partial2F

			ogsChem::LocalVector der_mob_R;
			der_mob_R  = LocalVector::Zero(_n_xi_Mob);
			NumDiff(_n_xi_Mob, delta_xi, std::bind(this->_ReductionGIA->Calc_Kin_Rate ,_1, loc_cur_xi_Sorp, loc_cur_xi_Sorp_tilde, loc_cur_xi_Sorp_bar, loc_cur_xi_Min, loc_cur_xi_Min_tilde,
					 loc_cur_xi_Min_bar, loc_cur_xi_Kin, loc_cur_xi_Kin_bar, loc_cur_eta, loc_cur_eta_bar, vec_Rate), vec_rate_old, loc_cur_xi_Mob, der_mob_R);

			ogsChem::LocalVector der_sorpB_R;
			der_sorpB_R  = LocalVector::Zero(_n_xi_Sorp_bar);
			NumDiff(_n_xi_Sorp, delta_xi, std::bind(this->_ReductionGIA->Calc_Kin_Rate ,loc_cur_xi_Mob, loc_cur_xi_Sorp, loc_cur_xi_Sorp_tilde, _4, loc_cur_xi_Min, loc_cur_xi_Min_tilde,
					 loc_cur_xi_Min_bar, loc_cur_xi_Kin, loc_cur_xi_Kin_bar, loc_cur_eta, loc_cur_eta_bar, vec_Rate), vec_rate_old, loc_cur_xi_Sorp_bar, der_sorpB_R);

			//ogsChem::LocalVector der_sorpB_li_R 	= der_sorpB_R.topRows(_n_xi_Sorp_bar_li);
			//ogsChem::LocalVector der_sorpB_ld_R 	= der_sorpB_R.bottomRows(_n_xi_Sorp_bar_ld);

			ogsChem::LocalVector der_minB_R;
			der_minB_R  = LocalVector::Zero(_n_xi_Min_bar);
			NumDiff(_n_xi_Sorp, delta_xi, std::bind(this->_ReductionGIA->Calc_Kin_Rate ,loc_cur_xi_Mob, loc_cur_xi_Sorp, loc_cur_xi_Sorp_tilde, loc_cur_xi_Sorp_bar, loc_cur_xi_Min, loc_cur_xi_Min_tilde,
					 _7, loc_cur_xi_Kin, loc_cur_xi_Kin_bar, loc_cur_eta, loc_cur_eta_bar, vec_Rate), vec_rate_old, loc_cur_xi_Min_bar, der_minB_R);

			mat_p2F.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde, 0, _n_xi_Sorp, _n_xi_Mob)							 = - dt * theta_water_content * _mat_Asorp * der_mob_R;
			mat_p2F.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp, 0, _n_xi_Min, _n_xi_Mob) 			 = - dt * theta_water_content * _mat_Amin * der_mob_R;
			mat_p2F.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp + _n_xi_Min, 0, _n_xi_Kin, _n_xi_Mob)  = - dt * theta_water_content * mat_A1kin * der_mob_R;

			//TODO IMPORTANT: li and ld is merged together. it should be checked later for its correctness.
			mat_p2F.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde, _n_xi_Mob, _n_xi_Sorp, _n_xi_Sorp_bar) 						 = - dt * theta_water_content * _mat_Asorp * der_sorpB_R;
			mat_p2F.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp, _n_xi_Mob, _n_xi_Min, _n_xi_Sorp_bar) 			 = - dt * theta_water_content * _mat_Amin * der_sorpB_R;
			mat_p2F.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp + _n_xi_Min, _n_xi_Mob, _n_xi_Kin, _n_xi_Sorp_bar) = - dt * theta_water_content * mat_A1kin * der_sorpB_R;

			mat_p2F.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde, _n_xi_Mob + _n_xi_Sorp_bar, _n_xi_Sorp, _n_xi_Min_bar) 						 = - dt * theta_water_content * _mat_Asorp * der_minB_R;
			mat_p2F.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp, _n_xi_Mob + _n_xi_Sorp_bar, _n_xi_Min, _n_xi_Min_bar) 		     = - dt * theta_water_content * _mat_Amin * der_minB_R;
			mat_p2F.block(_n_xi_Sorp_tilde + _n_xi_Min_tilde + _n_xi_Sorp + _n_xi_Min, _n_xi_Mob + _n_xi_Sorp_bar, _n_xi_Kin, _n_xi_Min_bar) = - dt * theta_water_content * mat_A1kin * der_minB_R;

			// add the relavent identity matrices
			mat_p2F.block(0, _n_xi_Mob, _n_xi_Sorp_bar_li, _n_xi_Sorp_bar_li)							 		 = LocalMatrix::Identity(_n_xi_Sorp_bar_li, _n_xi_Sorp_bar_li);
			mat_p2F.block(_n_xi_Sorp_bar_li, _n_xi_Mob + _n_xi_Sorp_bar_li, _n_xi_Min_bar, _n_xi_Sorp_bar_ld)	 = mat_Ald;
			mat_p2F.block(_n_xi_Sorp_bar_li, _n_xi_Mob + _n_xi_Sorp_bar, _n_xi_Min_bar, _n_xi_Min_bar)	 		 = LocalMatrix::Identity(_n_xi_Min_bar, _n_xi_Min_bar);

			////// calculate vprime
			Vprime(vec_conc, mat_vprime);

			// construct local Jacobian matrix
			Jacobian_local = mat_p1F + mat_p2F * mat_vprime;

			// construct global Jacobian matrix
			Jacobian_global.block(node_idx * _n_xi_global, node_idx * _n_xi_global, _n_xi_global, _n_xi_global) += Jacobian_local;

		}  // end of if statement
	}  // end of node based for loop






}



template <class T1, class T2>
void FunctionReductConc<T1, T2>::NumDiff(std::size_t & col,
			              			 	 ogsChem::LocalVector & delta_xi,
			              			 	 ogsChem::LocalVector & f,
			              			 	 ogsChem::LocalVector & f_old,
			              			 	 ogsChem::LocalVector & unknown,
			              			 	 ogsChem::LocalVector & DrateDxi)
{
	ogsChem::LocalVector xi = LocalVector::Zero(col);
for(std::size_t i = 0; i < col; i++)
{
    xi        = unknown;
    xi(i)     = xi(i) + delta_xi * xi(i).norm();
    DrateDxi(i) = ( f(xi) - f_old) / (delta_xi * xi(i).norm());  //is it elementwise?
}

}

template <class T1, class T2>
void FunctionReductConc<T1, T2>::Vprime( MathLib::LocalVector & vec_conc,
										 MathLib::LocalMatrix & mat_vprime)
{


}

