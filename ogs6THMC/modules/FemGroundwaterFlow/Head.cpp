
#include "Head.h"

#include "Ogs6FemData.h"


void FunctionHead::initialize(const BaseLib::Options &option)
{
    Ogs6FemData* femData = Ogs6FemData::getInstance();

    //const BaseLib::Options* op = option.getSubGroup("ProcessData")->getSubGroup("GROUNDWATER_FLOW");
    size_t msh_id = option.getOption<size_t>("MeshID");
    size_t time_id = option.getOption<size_t>("TimeGroupID");

    NumLib::ITimeStepFunction* tim = femData->list_tim[time_id];

    //mesh and FE objects
    DiscreteLib::DiscreteSystem* dis = femData->list_dis_sys[msh_id];
    MeshLib::IMesh* msh = dis->getMesh();
    _feObjects = new FemLib::LagrangianFeObjectContainer(*msh);

    // equations
    GWFemEquation::LinearAssemblerType* linear_assembler = new GWFemEquation::LinearAssemblerType(*_feObjects);
    GWFemEquation::ResidualAssemblerType* r_assembler = new GWFemEquation::ResidualAssemblerType(*_feObjects);
    GWFemEquation::JacobianAssemblerType* j_eqs = new GWFemEquation::JacobianAssemblerType(*_feObjects);
    GWFemEquation* eqs = new  GWFemEquation(linear_assembler, r_assembler, j_eqs);

    // set up problem
    _problem = new GWFemProblem(dis);
    _problem->setEquation(eqs);
    _problem->setTimeSteppingFunction(*tim);

    // set up solution
    _solution = new MySolutionType(dis, _problem);
    MySolutionType::LinearSolverType* linear_solver = _solution->getLinearEquationSolver();
    linear_solver->setOption(option);
    _solution->getNonlinearSolver()->setOption(option);

    this->setOutput(Head, _problem->getVariable(0)->getIC());
}
