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
#include "NestedLocalProbNRIterationStepInitializer.h"
#include "MathLib/ODE/RungeKutta4.h"

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
	// initialize xi_mob
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
	//MyNonLinearAssemblerType* non_linear_assembler = new MyNonLinearAssemblerType(_feObjects, this->_ReductionGIA );
	//MyNonLinearResidualAssemblerType* non_linear_r_assembler = new MyNonLinearResidualAssemblerType(_feObjects, this->_ReductionGIA );
	//MyNonLinearJacobianAssemblerType* non_linear_j_assembler = new MyNonLinearJacobianAssemblerType(_feObjects, this->_ReductionGIA, this);

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
//	_non_linear_problem = new MyNonLinearReactiveTransportProblemType(dis);
//	_non_linear_eqs     = _non_linear_problem->createEquation();
//	_non_linear_eqs->initialize( non_linear_assembler, non_linear_r_assembler, non_linear_j_assembler );
//	_non_linear_problem->setTimeSteppingFunction(*tim);
	// for nonlinear coupled transport problem, variables are xi_mobile species
//	for ( i=0; i < _n_xi_global ; i++ )
//	{
//		std::stringstream str_tmp;
//		str_tmp << "xi_global_" << i ;
//		_non_linear_problem->addVariable( str_tmp.str() );
//	}

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

	// set IC for eta_mob
	for ( i=0; i < _n_eta; i++ )
	{
		SolutionLib::FemIC* eta_ic = new SolutionLib::FemIC(msh);
		eta_ic->addDistribution( femData->geo->getDomainObj(), new NumLib::TXFunctionDirect<double>( _eta[i]->getDiscreteData() ) );
		_linear_problems[i]->getVariable(0)->setIC( eta_ic );  // be careful, here each linear problem only has one variable "eta".
	}
	// set IC for xi_mob
	for ( i=0; i < _n_xi_global; i++ )
	{
		SolutionLib::FemIC* xi_mob_ic = new SolutionLib::FemIC(msh);
		xi_mob_ic->addDistribution( femData->geo->getDomainObj(), new NumLib::TXFunctionDirect<double>( _xi_global[i]->getDiscreteData() ) );
		_non_linear_problem->getVariable(i)->setIC( xi_mob_ic );
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
//	myNRIterator = new MyNRIterationStepInitializer(non_linear_r_assembler, non_linear_j_assembler);
//	myNSolverFactory = new MyDiscreteNonlinearSolverFactory( myNRIterator );
//	this->_non_linear_solution = new MyNonLinearSolutionType( dis, this->_non_linear_problem, myNSolverFactory );
//    this->_non_linear_solution->getDofEquationIdTable()->setNumberingType(DiscreteLib::DofNumberingType::BY_POINT);  // global order
//    this->_non_linear_solution->getDofEquationIdTable()->setLocalNumberingType(DiscreteLib::DofNumberingType::BY_VARIABLE);  // local order
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

    // set initial output parameter
	for (i=0; i<_concentrations.size(); i++) {
		OutputVariableInfo var1(this->getOutputParameterName(i), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _concentrations[i]);
        femData->outController.setOutput(var1.name, var1);
    }

//#ifdef _DEBUG
//    // -----------debugging, output eta and xi----------------------
//    for (i=0; i < _n_eta_mob; i++) {
//        std::stringstream str_tmp;
//		str_tmp << "eta_mob_" << i ;
//        OutputVariableInfo var1(str_tmp.str(), _msh_id,  OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _eta_mob[i]);
//        femData->outController.setOutput(var1.name, var1);
//    }
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
//    for (i=0; i < _n_xi_immob; i++) {
//        std::stringstream str_tmp1, str_tmp2;
//		str_tmp1 << "xi_immob_" << i ;
//        OutputVariableInfo var1(str_tmp1.str(), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _xi_immob[i]);
//        femData->outController.setOutput(var1.name, var1);
//        str_tmp2 << "xi_immob_rate_" << i;
//        OutputVariableInfo var2(str_tmp2.str(), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _xi_immob_rates[i]);
//        femData->outController.setOutput(var2.name, var2);
//    }
//    // -----------------end of debugging-----------------------------
//#endif

    linear_solver = NULL;
    optNum = NULL;

    return true;
}

template <class T1, class T2>
void FunctionReductConc<T1, T2>::initializeTimeStep(const NumLib::TimeStep &/*time*/)
{
	//size_t i; 
 //   const NumLib::ITXFunction *vel = this->getInput<NumLib::ITXFunction>(Velocity);

	//// set velocity for linear problem
	//for ( i=0; i < _linear_problems.size(); i++ ) {
	//	_linear_problems[i]->getEquation()->getLinearAssembler()->setVelocity(vel);
	//	_linear_problems[i]->getEquation()->getResidualAssembler()->setVelocity(vel);
	//	_linear_problems[i]->getEquation()->getJacobianAssembler()->setVelocity(vel);
	//}
	//// set velocity for nonlinear problem as well
	//_non_linear_problem->getEquation()->getLinearAssembler()->setVelocity(vel); 
 //   _non_linear_problem->getEquation()->getResidualAssembler()->setVelocity(vel); 
	//_non_linear_problem->getEquation()->getJacobianAssembler()->setVelocity(vel); 
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
    //convert_eta_xi_to_conc(); 
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
//
//#ifdef _DEBUG
//    // -----------debugging, output eta and xi----------------------
//    for (i=0; i < _n_eta_mob; i++) {
//        std::stringstream str_tmp;
//		str_tmp << "eta_mob_" << i ;
//        OutputVariableInfo var1(str_tmp.str(), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _eta_mob[i]);
//        femData->outController.setOutput(var1.name, var1);
//    }
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
//#endif
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
	const double epsilon  = 1.0e-6;
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
	loc_eta         		  = LocalVector::Zero( _n_eta );
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
void FunctionReductConc<T1, T2>::calc_nodal_xi_immob_ode(double dt)
{
	size_t i, node_idx;

    double t0 = 0.0;

	// initialize the local vector
	MathLib::LocalVector loc_eta;
	MathLib::LocalVector loc_eta_bar;
	MathLib::LocalVector loc_xi_global;
	MathLib::LocalVector loc_xi_local;
    MathLib::LocalVector loc_xi_local_rate;
	MathLib::LocalVector loc_xi_local_new;

	// initialize the local vector
	loc_eta        		= LocalVector::Zero( _n_eta );
	loc_eta_bar         = LocalVector::Zero( _n_eta_bar );
	loc_xi_global       = LocalVector::Zero( _n_xi_global );
	loc_xi_local       	= LocalVector::Zero( _n_xi_local );
    loc_xi_local_rate  	= LocalVector::Zero( _n_xi_local );
	loc_xi_local_new   	= LocalVector::Zero( _n_xi_local );

	// signal, solving local ODEs of xi_immob
	INFO("--Solving local ODE problem of xi_immob...");

	// initialize the ODE Runge-Kutta solution class
	// MathLib::RungeKutta4<Local_ODE_Xi_immob, MathLib::LocalVector>* rk4 = new MathLib::RungeKutta4<Local_ODE_Xi_immob, MathLib::LocalVector>();
    // or using the StepperBulischStoer class
//    MathLib::StepperBulischStoer<Local_ODE_Xi_immob>* sbs
//        = new MathLib::StepperBulischStoer<Local_ODE_Xi_immob>(loc_xi_immob,
//                                                               loc_xi_immob_rate,
//                                                               t0,
//                                                               1.0e-6,
//                                                               1.0e-6,
//                                                               true);



//	// loop over all the nodes
//    for (node_idx = _concentrations[0]->getDiscreteData()->getRangeBegin();
//	     node_idx < _concentrations[0]->getDiscreteData()->getRangeEnd();
//		 node_idx++ )
//	{
//		// skip the boundary nodes
//		if ( ! this->_solution->isBCNode(node_idx) )
//		{
//			// on each node, get the right start value
//			// get the right set of eta and xi
//			for (i=0; i < _n_eta_mob; i++)
//				loc_eta_mob[i] = this->_eta_mob[i]->getValue(node_idx);
//			// fill in eta_immob
//			for (i=0; i < _n_eta_immob; i++)
//				loc_eta_immob[i] = this->_eta_immob[i]->getValue(node_idx);
//			for (i=0; i < _n_xi_mob; i++)
//				loc_xi_mob[i] = this->_xi_mob[i]->getValue(node_idx);
//			for (i=0; i < _n_xi_immob; i++)
//				loc_xi_immob[i] = this->_xi_immob[i]->getValue(node_idx); // HS, here are always the old xi_immob values.
//
//			// get the right reference values to ODE RHS function
//			this->_local_ode_xi_immob->update_eta_xi( loc_eta_mob, loc_eta_immob, loc_xi_mob, loc_xi_immob);
//            loc_xi_immob_rate = (*_local_ode_xi_immob)(dt, loc_xi_immob);
//            sbs->set_y(loc_xi_immob);
//            sbs->set_dydx(loc_xi_immob_rate);
//
//			// solve the local ODE problem using RK
//			// rk4->solve( *_local_ode_xi_immob, 0.0, dt, loc_xi_immob, loc_xi_immob_new);
//            sbs->step( dt, _local_ode_xi_immob );
//            loc_xi_immob_new = sbs->get_y();
//
//			// collect the xi_immob_new
//			for (i=0; i < _n_xi_immob; i++)
//				_xi_immob_new[i]->setValue(node_idx, loc_xi_immob_new[i]);
//		} // end of if
//        else
//        {
//            // if it is boundary nodes, then xi_immob_new is equal to xi_immob.
//            for (i=0; i < _n_xi_immob; i++)
//				_xi_immob_new[i]->setValue(node_idx, _xi_immob[i]->getValue( node_idx ) );
//        }
//
//	}  // end of for node_idx

    // delete rk4;
//    delete sbs;
}
