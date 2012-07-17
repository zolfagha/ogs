
#include "Head.h"

#include "Ogs6FemData.h"

//namespace Geo
//{

void FunctionHead::initialize(const BaseLib::Options &option)
{
	const BaseLib::Options* op = option.getSubGroup("ProcessData")->getSubGroup("GROUNDWATER_FLOW");
	std::string msh_key = op->getOption("Mesh");
	std::string mat_key = op->getOption("Material");

	DiscreteLib::DiscreteSystem* dis = 0;
	MaterialLib::PorousMedia *pm = 0;

    _feObjects = new FemLib::LagrangianFeObjectContainer(*dis->getMesh());
    GWFemEquation::LinearAssemblerType* linear_assembler = new GWFemEquation::LinearAssemblerType(*_feObjects, *pm);
    GWFemEquation::ResidualAssemblerType* r_assembler = new GWFemEquation::ResidualAssemblerType(*_feObjects, *pm);
    GWFemEquation::JacobianAssemblerType* j_eqs = new GWFemEquation::JacobianAssemblerType(*_feObjects, *pm);
    GWFemEquation* eqs = new  GWFemEquation(linear_assembler, r_assembler, j_eqs);
    GWFemProblem* _problem = new GWFemProblem(dis);
    _problem->setEquation(eqs);
	_solution = new MySolutionType(dis, _problem);
    //_solHead->getTimeODEAssembler()->setTheta(1.0);
	MySolutionType::LinearSolverType* linear_solver = _solution->getLinearEquationSolver();
    linear_solver->setOption(option);
    _solution->getNonlinearSolver()->setOption(option);

    this->setOutput(Head, _problem->getVariable(0)->getIC());
}

//OGS_PROCESS(GROUNDWATER_FLOW, FunctionHead);

//}
