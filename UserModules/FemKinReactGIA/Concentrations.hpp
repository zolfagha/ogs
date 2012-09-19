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

#include "DiscreteLib/Utils/DiscreteSystemContainerPerMesh.h"
#include "FemLib/Function/FemNodalFunction.h"
#include "NumLib/Function/TXFunctionBuilder.h"
#include "OutputIO/OutputBuilder.h"
#include "OutputIO/OutputTimingBuilder.h"
#include "SolutionLib/Fem/FemSourceTerm.h"
#include "Ogs6FemData.h"

template <class T1, class T2>
bool FunctionConcentrations<T1,T2>::initialize(const BaseLib::Options &option)
{
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    size_t msh_id = option.getOption<size_t>("MeshID");
    size_t time_id = option.getOption<size_t>("TimeGroupID");
    NumLib::ITimeStepFunction* tim = femData->list_tim[time_id];

    //mesh and FE objects
    MeshLib::IMesh* msh = femData->list_mesh[msh_id];
    MyDiscreteSystem* dis = 0;
    dis = DiscreteLib::DiscreteSystemContainerPerMesh::getInstance()->createObject<MyDiscreteSystem>(msh);
    _feObjects = new FemLib::LagrangianFeObjectContainer(*msh);

	// first get the number of components
	size_t n_Comp = femData->map_ChemComp.size(); 
	// add all concentrations to discretized memory space
	NumLib::LocalVector local_conc(n_Comp);
    local_conc *= .0;
	_concentrations = new MyNodalFunctionVector(); 
	_concentrations->initialize(*dis, FemLib::PolynomialOrder::Linear, local_conc);
	
	// get the transformation class instance here
	this->_ReductionKin = femData->m_KinReductScheme; 
	// make sure the reduction scheme is already initialized. 
	if ( !(this->_ReductionKin->IsInitialized()) ) 
	{
		// error msg
	    ERR("While initialize the Global Implicit Reactive Transport Process, the reduction scheme has not been correctly initialized! ");
		// then stop the program
		exit(1);
	}

	// tell me how many eta and how many xi we have
	size_t n_eta, n_xi, n_eta_mob; 
	// get n_eta and n_xi
	n_eta     = this->_ReductionKin->get_n_eta();
	n_eta_mob = this->_ReductionKin->get_n_eta_mob(); 
	n_xi      = this->_ReductionKin->get_n_xi(); 
	// based on the transformation class instance, add eta and xi to the discretized memory space
	NumLib::LocalVector local_eta(n_eta), local_xi(n_xi);
	local_eta *= .0; local_xi *= .0; 
	this->_eta = new MyNodalFunctionVector(); 
	this->_eta ->initialize( *dis, FemLib::PolynomialOrder::Linear, local_eta );
    this->_xi  = new MyNodalFunctionVector(); 
	this->_xi  ->initialize( *dis, FemLib::PolynomialOrder::Linear, local_xi  );

	// eta goes to the linear problems

	// xi  goes to the nonlinear problems

	// first set up the linear problem for eta
	// equations
    MyLinearAssemblerType* linear_assembler = new MyLinearAssemblerType(_feObjects);
    MyResidualAssemblerType* r_assembler = new MyResidualAssemblerType(_feObjects);
    MyJacobianAssemblerType* j_eqs = new MyJacobianAssemblerType(_feObjects);

	// set up problem

	// set up variables

	// IC

	// BC

	// ST

	// set up solution

	// set initial output

	// initial output parameter



	/*
    // equations
    MyLinearAssemblerType* linear_assembler = new MyLinearAssemblerType(_compound, _feObjects);
    MyResidualAssemblerType* r_assembler = new MyResidualAssemblerType(_compound, _feObjects);
    MyJacobianAssemblerType* j_eqs = new MyJacobianAssemblerType(_compound, _feObjects);

    // set up problem
    _problem = new MyProblemType(dis);
    MyEquationType* eqs = _problem->createEquation();
    eqs->initialize(linear_assembler, r_assembler, j_eqs);
    _problem->setTimeSteppingFunction(*tim);

    // set up variables
	// in this case, the variables includes: \
	// 1) concentrations of all components in the MCP data structure
    // typename MyProblemType::MyVariable* concentrations = _problem->addVariable("concentrations");
	for ( size_t i=0; i < n_Comp ; i++ )
		_problem->addVariable( femData->map_ChemComp[i]->second->get_name() );
	// 2) values of transformed eta

	// 3) values of transformed xi


    // IC
    NumLib::TXFunctionBuilder f_builder;
    typename MyProblemType::MyVariable::MyNodalFunctionScalar* c0 = new typename MyProblemType::MyVariable::MyNodalFunctionScalar();
    c0->initialize(*dis, FemLib::PolynomialOrder::Linear, 0);
    concentrations->setIC(c0);

    // BC
    const BaseLib::Options* opBCList = option.getSubGroup("BCList");
    for (const BaseLib::Options* opBC = opBCList->getFirstSubGroup("BC"); opBC!=0; opBC = opBCList->getNextSubGroup())
    {
        std::string geo_type = opBC->getOption("GeometryType");
        std::string geo_name = opBC->getOption("GeometryName");
        const GeoLib::GeoObject* geo_obj = femData->geo->searchGeoByName(femData->geo_unique_name, geo_type, geo_name);
        std::string dis_name = opBC->getOption("DistributionType");
        double dis_v = opBC->getOption<double>("DistributionValue");
        NumLib::ITXFunction* f_bc =  f_builder.create(dis_name, dis_v);
        concentrations->addDirichletBC(new SolutionLib::FemDirichletBC(msh, geo_obj, f_bc));
    }

    // ST
    const BaseLib::Options* opSTList = option.getSubGroup("STList");
    for (const BaseLib::Options* opST = opSTList->getFirstSubGroup("ST"); opST!=0; opST = opSTList->getNextSubGroup())
    {
        std::string geo_type = opST->getOption("GeometryType");
        std::string geo_name = opST->getOption("GeometryName");
        const GeoLib::GeoObject* geo_obj = femData->geo->searchGeoByName(femData->geo_unique_name, geo_type, geo_name);
        std::string st_type = opST->getOption("STType");
        std::string dis_name = opST->getOption("DistributionType");
        double dis_v = opST->getOption<double>("DistributionValue");
        if (st_type.compare("NEUMANN")==0) {
            dis_v *= -1; // user set inflow as positive sign but internally negative
        }
        NumLib::ITXFunction* f_st =  f_builder.create(dis_name, dis_v);
        if (f_st!=NULL) {
            SolutionLib::IFemNeumannBC *femSt = 0;
            if (st_type.compare("NEUMANN")==0) {
                femSt = new SolutionLib::FemNeumannBC(msh, _feObjects, geo_obj, f_st);
            } else if (st_type.compare("SOURCESINK")==0) {
                femSt = new SolutionLib::FemSourceTerm(msh, geo_obj, f_st);
            }
            concentrations->addNeumannBC(femSt);
        } else {
            WARN("Distribution type %s is specified but not found. Ignore this ST.", dis_name.c_str());
        }
    }

    // set up solution
    _solution = new MySolutionType(dis, _problem);
    MyLinearSolver* linear_solver = _solution->getLinearEquationSolver();
    const BaseLib::Options* optNum = option.getSubGroup("Numerics");
    linear_solver->setOption(*optNum);
    _solution->getNonlinearSolver()->setOption(*optNum);

    // set initial output
    OutputVariableInfo var(this->getOutputParameterName(Concentrations), OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _solution->getCurrentSolution(0));
    femData->outController.setOutput(var.name, var); 

    // initial output parameter
    this->setOutput(Concentrations, concentrations->getIC());
	*/

    return true;
}

template <class T1, class T2>
void FunctionConcentrations<T1, T2>::initializeTimeStep(const NumLib::TimeStep &/*time*/)
{
    const NumLib::ITXFunction *vel = this->getInput<NumLib::ITXFunction>(Velocity);
	// set velocity for linear problem
	this->_linear_problem->getEquation()->getLinearAssembler()->setVelocity(vel);
    this->_linear_problem->getEquation()->getResidualAssembler()->setVelocity(vel);
    this->_linear_problem->getEquation()->getJacobianAssembler()->setVelocity(vel);
	// set velocity for nonlinear problem as well
	// TODO
}

template <class T1, class T2>
void FunctionConcentrations<T1, T2>::updateOutputParameter(const NumLib::TimeStep &/*time*/)
{
	// TODO 
    // setOutput(Concentrations, _solution->getCurrentSolution(0));
}

template <class T1, class T2>
void FunctionConcentrations<T1, T2>::output(const NumLib::TimeStep &/*time*/)
{
    //update data for output
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    
	// TODO the following needs to be changed. 
	// OutputVariableInfo var(this->getOutputParameterName(Concentrations), OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _solution->getCurrentSolution(0));
    // femData->outController.setOutput(var.name, var); 
};
