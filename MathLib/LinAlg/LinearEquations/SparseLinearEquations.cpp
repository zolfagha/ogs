
#include "SparseLinearEquations.h"

#include <iostream>
#include <algorithm>

#include "MathLib/LinAlg/Solvers/CG.h"
#include "MathLib/LinAlg/Solvers/BiCGStab.h"
#include "MathLib/sparse.h"

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

void SparseLinearEquations::solveEqs(CRSMatrix<double, unsigned> *A, double *rhs, double *x)
{
    double eps = _option.error_tolerance;
    size_t steps =  _option.max_iteration_step;
    switch (_option.solver_type)
    {
    case SolverCG:
        CG(A, rhs, x, eps, steps);
        std::cout << "MathLib::CG converged within " << steps << ", residuum is " << eps << std::endl;
        break;
    case SolverBiCGStab:
        BiCGStab(*A, rhs, x, eps, steps);
    default:
        break;
    }
}

} // end namespace

