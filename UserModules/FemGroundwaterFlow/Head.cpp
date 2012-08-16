/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Head.cpp
 *
 * Created on 2012-07-13 by Norihiro Watanabe
 */

#include "Head.h"

#include "FemLib/Function/FemNodalFunction.h"
#include "NumLib/Function/TXFunctionBuilder.h"
#include "OutputIO/OutputBuilder.h"
#include "OutputIO/OutputTimingBuilder.h"
#include "SolutionLib/Fem/FemSourceTerm.h"
#include "Ogs6FemData.h"

bool FunctionHead::initialize(const BaseLib::Options &option)
{
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    size_t msh_id = option.getOption<size_t>("MeshID");
    size_t time_id = option.getOption<size_t>("TimeGroupID");
    NumLib::ITimeStepFunction* tim = femData->list_tim[time_id];

    //mesh and FE objects
    MyProblemType::MyDiscreteSystem* dis = (MyProblemType::MyDiscreteSystem*)femData->list_dis_sys[msh_id]; //TODO
    MeshLib::IMesh* msh = dis->getMesh();
    _feObjects = new FemLib::LagrangianFeObjectContainer(*msh);

    // local assemblers
    FemGWEquation::LinearAssemblerType* linear_assembler = new FemGWEquation::LinearAssemblerType(*_feObjects);
    FemGWEquation::ResidualAssemblerType* r_assembler = new FemGWEquation::ResidualAssemblerType(*_feObjects);
    FemGWEquation::JacobianAssemblerType* j_eqs = new FemGWEquation::JacobianAssemblerType(*_feObjects);

    // set up problem
    _problem = new FemGWProblem(dis);
    FemGWEquation* eqs = _problem->createEquation();
    eqs->initialize(linear_assembler, r_assembler, j_eqs);
    _problem->setTimeSteppingFunction(*tim);
    // set up variable
    MyProblemType::MyVariable* head = _problem->addVariable("head");
    // IC
    NumLib::TXFunctionBuilder f_builder;
    MyProblemType::MyVariable::MyNodalFunctionScalar* h0 = new MyProblemType::MyVariable::MyNodalFunctionScalar();
    h0->initialize(*dis, FemLib::PolynomialOrder::Linear, 0);
    head->setIC(h0);
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
        head->addDirichletBC(new SolutionLib::FemDirichletBC(msh, geo_obj, f_bc));
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
            head->addNeumannBC(femSt);
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
    OutputVariableInfo var("HEAD", OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _solution->getCurrentSolution(0));
    femData->outController.setOutput("HEAD", var); 

    // initial output parameter
    this->setOutput(Head, head->getIC());

    return true;
}

void FunctionHead::updateOutputParameter(const NumLib::TimeStep &time)
{
    setOutput(Head, _solution->getCurrentSolution(0));
}

void FunctionHead::output(const NumLib::TimeStep &time) 
{
    //update data for output
    Ogs6FemData* femData = Ogs6FemData::getInstance();
    OutputVariableInfo var("HEAD", OutputVariableInfo::Node, OutputVariableInfo::Real, 1, _solution->getCurrentSolution(0));
    femData->outController.setOutput("HEAD", var); 
};
