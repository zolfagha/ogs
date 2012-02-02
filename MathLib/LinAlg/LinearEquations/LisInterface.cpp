
#include "LisInterface.h"

#include <iostream>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "lis.h"

namespace MathLib
{

void LIS_Solver::initialize(int argc, char *argv[])
{
    lis_initialize(&argc, &argv);
}

void LIS_Solver::finialize()
{
    lis_finalize();
}

void LIS_Solver::solve(CRSSigned *A, double *x, double *b, LIS_option &option)
{
    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "*** LIS solver computation" << std::endl;

    // Creating a matrix.
    LIS_MATRIX AA;
    LIS_VECTOR bb,xx;
    LIS_SOLVER solver;
    int ierr = lis_matrix_create(0, &AA);
    ierr = lis_matrix_set_type(AA, LIS_MATRIX_CRS);
    ierr = lis_matrix_set_size(AA, 0, A->dimension);

    // Matrix solver and Precondition can be handled better way.
    const size_t MAX_ZEILE = 512;
    char solver_options[MAX_ZEILE], tol_option[MAX_ZEILE];

    int nthreads = omp_get_num_threads();
    //omp_set_num_threads (nthreads);

    sprintf(solver_options, "-i %d -p %d %s", option.ls_method, option.ls_precond, option.ls_extra_arg.c_str()); 
    sprintf(tol_option, "-tol %e -maxiter %d -omp_num_threads %d -initx_zeros 0", option.ls_error_tolerance, option.ls_max_iterations, nthreads);

    ierr = lis_matrix_set_crs(A->nonzero, A->row_ptr, A->col_idx, A->data, AA);
    ierr = lis_matrix_assemble(AA);

    // Assemble the vector, b, x
    ierr = lis_vector_duplicate(AA, &bb);
    ierr = lis_vector_duplicate(AA, &xx);
#pragma omp parallel for
    for (long i=0; i < A->dimension; ++i)
    {
        ierr = lis_vector_set_value(LIS_INS_VALUE, i, x[i], xx);
        ierr = lis_vector_set_value(LIS_INS_VALUE, i, b[i], bb);
    }

    // Create solver
    ierr = lis_solver_create(&solver);

    ierr = lis_solver_set_option(solver_options, solver);
    ierr = lis_solver_set_option(tol_option, solver);
    ierr = lis_solver_set_option("-print mem", solver);
    
    ierr = lis_solve(AA, bb, xx, solver);
    int iter = 0;
    ierr = lis_solver_get_iters(solver, &iter);
    //NW
    printf("\t iteration: %d/%d\n", iter, option.ls_max_iterations);
    double resid = 0.0;
    ierr = lis_solver_get_residualnorm(solver, &resid);
    printf("\t residuals: %e\n", resid);
    //	lis_vector_print(xx);
    //	lis_vector_print(bb);

    // Update the solution (answer) into the x vector
#pragma omp parallel for
    for(long i=0; i<A->dimension; ++i)
    {
        lis_vector_get_value(xx,i,&(x[i]));
    }

    // Clear memory
    //lis_matrix_destroy(AA);
    lis_vector_destroy(bb);
    lis_vector_destroy(xx);
    lis_solver_destroy(solver);
    std::cout << "------------------------------------------------------------------" << std::endl;
}


}
