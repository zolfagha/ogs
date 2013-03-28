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
#include "MathLib/ODE/RungeKutta4.h"

template <class T1, class T2>
bool FunctionOPSConc<T1,T2>::initialize(const BaseLib::Options &option)
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
    this->_local_eq_react_sys = femData->m_EqReactSys; 
    // make sure the reduction scheme is already initialized. 
    if ( !(this->_local_eq_react_sys->IsInitialized()) ) 
	{
		// error msg
	    ERR("While initialize the OPS Reactive Transport Process, the chemEqReactSys class has not been correctly initialized! ");
		// then stop the program
		exit(1);
	}

    // first get the number of components
    _n_Comp = this->_local_eq_react_sys->get_n_Comp(); 

	// set concentrations of all components as output
	for ( i=0; i < _n_Comp; i++ )
		this->setOutputParameterName( i, femData->map_ChemComp[i]->second->get_name() ); 

	// creating local memory space to store IC and BC
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

	// linear assemblers
    MyLinearAssemblerType* linear_assembler = new MyLinearAssemblerType(_feObjects);
    MyLinearResidualAssemblerType* linear_r_assembler = new MyLinearResidualAssemblerType(_feObjects);
    MyLinearJacobianAssemblerType* linear_j_eqs = new MyLinearJacobianAssemblerType(_feObjects);


//	// for the linear transport problem, variables are eta_mobile
//	for ( i=0; i < _n_eta_mob ; i++ )
//	{
//		// set up problem
//		MyLinearTransportProblemType* linear_problem = new MyLinearTransportProblemType(dis);
//		MyLinearEquationType* linear_eqs = linear_problem->createEquation();
//		linear_eqs->initialize(linear_assembler, linear_r_assembler, linear_j_eqs);
//		linear_problem->setTimeSteppingFunction(*tim);
//		
//		// set up variables
//		// in this case, the variables includes: eta_0, eta_1, eta_2......, 
//		// which are the concentrations of all eta_mobile 
//		std::stringstream str_tmp;
//		str_tmp << "eta_" << i ; 
//		linear_problem->addVariable( str_tmp.str() );
//		_linear_problems.push_back(linear_problem); 
//	}
//	

//	// reduction problem
//	_problem = new MyKinReductionProblemType( dis, _ReductionKin ); 
//	_problem->setTimeSteppingFunction(*tim);  // applying the same time stepping function for all linear, non-linear and reduction problems
//	// add variables to the KinReduction problem class
//	// for the KinReduction problem, variables are the concentrations of all chemical components
//	// add all concentrations to discretized memory space
//	for ( i=0; i < _n_Comp; i++ )
//	{
//		MyVariableConc* comp_conc = _problem->addVariable( femData->map_ChemComp[i]->second->get_name() );
//        FemVariableBuilder var_builder;
//        var_builder.doit(femData->map_ChemComp[i]->second->get_name(), option, msh, femData->geo, femData->geo_unique_name, _feObjects, comp_conc);
//	}


    // set up linear solution
	for ( i=0; i < _n_Comp; i++ )
	{
		MyLinearSolutionType* linear_solution = new MyLinearSolutionType( dis, this->_linear_problems[i] ); 
		MyLinearSolver* linear_solver = linear_solution->getLinearEquationSolver();
		const BaseLib::Options* optNum = option.getSubGroup("Numerics");
		linear_solver->setOption(*optNum);
		this->_linear_solutions.push_back( linear_solution ); 
	}


	// set up solution
    // _solution = new MyKinReductionSolution(dis, _problem, this, _linear_problems, _linear_solutions, _non_linear_problem, _non_linear_solution);
//	
//    this->setOutput(Concentrations, _solution->getCurrentSolution(0));
//
//    // set initial output parameter
//	for (i=0; i<_concentrations.size(); i++) {
//		OutputVariableInfo var1(this->getOutputParameterName(i), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _concentrations[i]);
//        femData->outController.setOutput(var1.name, var1);
//    }
//
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
//
//    linear_solver = NULL;
//    optNum = NULL;

    return true;
}

template <class T1, class T2>
void FunctionOPSConc<T1, T2>::initializeTimeStep(const NumLib::TimeStep &/*time*/)
{
	size_t i; 
    const NumLib::ITXFunction *vel = this->getInput<NumLib::ITXFunction>(Velocity);

	// set velocity for linear problem
	for ( i=0; i < _linear_problems.size(); i++ ) {
		_linear_problems[i]->getEquation()->getLinearAssembler()->setVelocity(vel);
		_linear_problems[i]->getEquation()->getResidualAssembler()->setVelocity(vel);
		_linear_problems[i]->getEquation()->getJacobianAssembler()->setVelocity(vel);
	}
}

template <class T1, class T2>
void FunctionOPSConc<T1, T2>::updateOutputParameter(const NumLib::TimeStep &/*time*/)
{ 
    // convert eta and xi back to concentrations
    //convert_eta_xi_to_conc(); 
}

template <class T1, class T2>
void FunctionOPSConc<T1, T2>::output(const NumLib::TimeStep &/*time*/)
{
//    size_t i;
//
//    // update data for output
//    Ogs6FemData* femData = Ogs6FemData::getInstance();
//    this->_msh_id = this->_problem->getDiscreteSystem()->getMesh()->getID(); 
//	// set the new output
//	for (i=0; i<_concentrations.size(); i++) {
//		OutputVariableInfo var1(this->getOutputParameterName(i), _msh_id, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _concentrations[i]);
//        femData->outController.setOutput(var1.name, var1);
//    }
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

