
#include "SparseLinearEquations.h"

#include <iostream>

#include "MathLib/LinAlg/Solvers/CG.h"
#include "MathLib/LinAlg/Solvers/BiCGStab.h"


namespace MathLib
{

void SparseLinearEquations::setOption(const Base::Options &option)
{
    const Base::Options *op = option.getSubGroup("SpLinearOptions");
    if (op==0) return;

    if (op->hasOption("solver_type"))
        _option.solver_type = (SolverType)op->getOptionAsNum<int>("solver_type");
    if (op->hasOption("precon_type"))
        _option.precon_type = (PreconditionerType)op->getOptionAsNum<int>("precon_type");
    if (op->hasOption("error_tolerance"))
        _option.error_tolerance = op->getOptionAsNum<double>("error_tolerance");
    if (op->hasOption("max_iteration_step"))
        _option.max_iteration_step = op->getOptionAsNum<int>("max_iteration_step");
}

void SparseLinearEquations::solve()
{
    double eps = _option.error_tolerance;
    size_t steps =  _option.max_iteration_step;

    switch (_option.solver_type)
    {
    case SolverCG:
        CG(getA(), getRHS(), getX(), eps, steps);
        std::cout << "MathLib::CG converged within " << steps << ", residuum is " << eps << std::endl;
        break;
    case SolverBiCGStab:
        BiCGStab(*getA(), getRHS(), getX(), eps, steps);
    default:
        break;
    }
}


} // end namespace

