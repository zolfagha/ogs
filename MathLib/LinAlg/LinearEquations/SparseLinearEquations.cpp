
#include "SparseLinearEquations.h"

#include <iostream>

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

void SparseLinearEquations::solve()
{
    if (_vec_knownX_id.size()>0) {
        CRSMatrix<double, unsigned>* tmp_A = getA()->clone();
        double *org_eqsRHS = getRHS();
        double *org_eqsX = getX();
        std::vector<double> _tmp_b;
        std::vector<double> _tmp_x;
        std::map<size_t,size_t> _map_solved_orgEqs;

        //std::cout << "#before\n";
        //_tmp_A->printMat();
        setKnownXi_ReduceSizeOfEQS(tmp_A, org_eqsRHS, org_eqsX, _vec_knownX_id, _vec_knownX_x, _tmp_b, _tmp_x, _map_solved_orgEqs);
        //std::cout << "\n#after\n";
        //_tmp_A->printMat();

        solveEqs(tmp_A, &_tmp_b[0], &_tmp_x[0], _option);

        const size_t dim = tmp_A->getNRows();
        for (size_t i=0; i<dim; i++) {
            setX(_map_solved_orgEqs[i], _tmp_x[i]);
        }

        delete tmp_A;
    } else {
        solveEqs(getA(), getRHS(), getX(), _option);
    }

}

void SparseLinearEquations::solveEqs(CRSMatrix<double, unsigned> *A, double *rhs, double *x, SpLinearOptions &option)
{
    double eps = option.error_tolerance;
    size_t steps =  option.max_iteration_step;
    switch (option.solver_type)
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

void SparseLinearEquations::setKnownXi_ReduceSizeOfEQS(CRSMatrix<double, unsigned> *A, double *org_eqsRHS, double *org_eqsX, const std::vector<size_t> &vec_id, const std::vector<double> &vec_x, std::vector<double> &out_b, std::vector<double> &out_x, std::map<size_t,size_t> &map_solved_orgEqs)
{
    assert(vec_id.size()==vec_x.size());

    const size_t n_org_rows = A->getNRows();

    std::vector<size_t> removed_rows(vec_id.size());
    for (size_t i=0; i<vec_id.size(); i++) {
        const size_t id = vec_id[i];
        const double val = vec_x[i];
        removed_rows[i] = id;

        //b_i -= A(i,k)*val, i!=k
        for (size_t j=0; j<A->getNCols(); j++)
            org_eqsRHS[j] -= A->getValue(j, id)*val;
        //b_k = A_kk*val
        org_eqsRHS[id] = val; //=eqsA(id, id)*val;
        org_eqsX[id] = val; //=eqsA(id, id)*val;
    }

    //remove rows and columns
    std::sort(removed_rows.begin(), removed_rows.end());
    A->eraseEntries(removed_rows.size(), &removed_rows[0]);
    const size_t n_new_rows = n_org_rows-removed_rows.size();

    //remove X,RHS
    out_b.resize(n_new_rows);
    out_x.resize(n_new_rows);
    size_t new_id = 0;
    for (size_t i=0; i<n_org_rows; i++) {
        if (std::find(removed_rows.begin(), removed_rows.end(), i)!=removed_rows.end()) continue;
        out_b[new_id] = org_eqsRHS[i];
        out_x[new_id] = org_eqsX[i];
        map_solved_orgEqs[new_id] = i;
        new_id++;
    }
}

} // end namespace

