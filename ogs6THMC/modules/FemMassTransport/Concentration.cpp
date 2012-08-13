/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Concentration.cpp
 *
 * Created on 2012-07-13 by Norihiro Watanabe
 */


#include "Concentration.h"

#include "logog.hpp"

#include "FemLib/Function/FemNodalFunction.h"
#include "NumLib/Function/TXFunctionBuilder.h"
#include "OutputIO/OutputBuilder.h"
#include "OutputIO/OutputTimingBuilder.h"
#include "Ogs6FemData.h"

bool FunctionConcentration::initialize(const BaseLib::Options &option)
{
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    size_t msh_id = option.getOption<size_t>("MeshID");
    size_t time_id = option.getOption<size_t>("TimeGroupID");
    NumLib::ITimeStepFunction* tim = femData->list_tim[time_id];

    //mesh and FE objects
    DiscreteLib::DiscreteSystem* dis = femData->list_dis_sys[msh_id];
    MeshLib::IMesh* msh = dis->getMesh();
    _feObjects = new FemLib::LagrangianFeObjectContainer(*msh);

    //Compound
    size_t compound_id = option.getOption<size_t>("CompoundID");
    if (femData->list_compound.size() < compound_id) {
        ERR("***Error: compound data not found (compound id=%d)", compound_id);
        return false;
    }
    _compound = femData->list_compound[compound_id];

    // equations
    MyEquationType::LinearAssemblerType* linear_assembler = new MyEquationType::LinearAssemblerType(_compound, _feObjects);
    MyEquationType::ResidualAssemblerType* r_assembler = new MyEquationType::ResidualAssemblerType(_compound, _feObjects);
    MyEquationType::JacobianAssemblerType* j_eqs = new MyEquationType::JacobianAssemblerType(_compound, _feObjects);
    MyEquationType* eqs = new  MyEquationType(linear_assembler, r_assembler, j_eqs);

    // set up problem
    _problem = new MyProblemType(dis);
    _problem->setEquation(eqs);
    _problem->setTimeSteppingFunction(*tim);
    // set up variable
    SolutionLib::FemVariable* concentration = _problem->addVariable("concentration");
    // IC
    NumLib::TXFunctionBuilder f_builder;
    FemLib::FemNodalFunctionScalar* c0 = new FemLib::FemNodalFunctionScalar(*dis, FemLib::PolynomialOrder::Linear, 0);
    concentration->setIC(c0);
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
        concentration->addDirichletBC(new SolutionLib::FemDirichletBC(msh, geo_obj, f_bc));
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
            concentration->addNeumannBC(femSt);
        } else {
            WARN("Distribution type %s is specified but not found. Ignore this ST.", dis_name.c_str());
        }
    }

    // set up solution
    _solution = new MySolutionType(dis, _problem);
    MySolutionType::LinearSolverType* linear_solver = _solution->getLinearEquationSolver();
    const BaseLib::Options* optNum = option.getSubGroup("Numerics");
    linear_solver->setOption(*optNum);
    _solution->getNonlinearSolver()->setOption(*optNum);

    // set initial output
    OutputVariableInfo var(_compound->name, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _solution->getCurrentSolution(0));
    femData->outController.setOutput(var.name, var); 

    // initial output parameter
    this->setOutput(Concentration, concentration->getIC());

    return true;
}

void FunctionConcentration::initializeTimeStep(const NumLib::TimeStep &/*time*/)
{
    const NumLib::ITXFunction *vel = this->getInput<NumLib::ITXFunction>(Velocity);
    this->_problem->getEquation()->getLinearAssembler()->setVelocity(vel);
    this->_problem->getEquation()->getResidualAssembler()->setVelocity(vel);
    this->_problem->getEquation()->getJacobianAssembler()->setVelocity(vel);
}

void FunctionConcentration::updateOutputParameter(const NumLib::TimeStep &/*time*/)
{
    setOutput(Concentration, _solution->getCurrentSolution(0));
}

void FunctionConcentration::output(const NumLib::TimeStep &/*time*/)
{
    //update data for output
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    OutputVariableInfo var(_compound->name, OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _solution->getCurrentSolution(0));
    femData->outController.setOutput(var.name, var); 
};
