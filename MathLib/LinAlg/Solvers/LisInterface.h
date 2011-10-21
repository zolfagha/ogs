
#pragma once

#include <string>
//#ifdef LIS
#include "lis.h"
//#endif
#include "LinAlg/Sparse/CRSMatrix.h"

/*
LIS matrix solver options
CG -i {cg|1}
BiCG -i {bicg|2}
CGS -i {cgs|3}
BiCGSTAB -i {bicgstab|4}
BiCGSTAB(l) -i {bicgstabl|5} -ell [2] Value for l
GPBiCG -i {gpbicg|6}
TFQMR -i {tfqmr|7}
Orthomin(m) -i {orthomin|8} -restart [40] Value for Restart m
GMRES(m) -i {gmres|9} -restart [40] Value for Restart m
Jacobi -i {jacobi|10}
Gauss-Seidel -i {gs|11}
SOR -i {sor|12} -omega [1.9] Value for Relaxation Coefficient  (0 <  < 2)
BiCGSafe -i {bicgsafe|13}
CR -i {cr|14}
BiCR -i {bicr|15}
CRS -i {crs|16}
BiCRSTAB -i {bicrstab|17}
GPBiCR -i {gpbicr|18}
BiCRSafe -i {bicrsafe|19}
FGMRES(m) -i {fgmres|20} -restart [40] Value for Restart m
IDR(s) -i {idrs|21} -restart [40] Value for Restart s

Preconditioner Option Auxiliary Option
None -p {none|0}
Jacobi -p {jacobi|1}
ILU(k) -p {ilu|2} -ilu_fill [0] Fill level k
SSOR -p {ssor|3} -ssor_w [1.0] Relaxation Coefficient  (0 <  < 2)
Hybrid -p {hybrid|4} -hybrid_i [sor] Iterative method
-hybrid_maxiter [25] Maximum number of iterations
-hybrid_tol [1.0e-3] Convergence criteria
-hybrid_w [1.5] Relaxation Coefficient  for
the SOR method (0 <  < 2)
-hybrid_ell [2] Value for l of the BiCGSTAB(l) method
-hybrid_restart [40] Restart values for GMRES and Orthomin
I+S -p {is|5} -is_alpha [1.0] Parameter ?for preconditioner
of a I + ?(m) type
-is_m [3] Parameter m for preconditioner
of a I + ?(m) type
SAINV -p {sainv|6} -sainv_drop [0.05] Drop criteria
SA-AMG -p {saamg|7} -saamg_unsym [false] Selection of asymmetric version
Crout ILU -p {iluc|8} -iluc_drop [0.05] Drop criteria
-iluc_rate [5.0] Ratio of Maximum fill-in
ILUT -p {ilut|9} -ilut_drop [0.05] Drop criteria
-ilut_rate [5.0] Ratio of Maximum fill-in
additive Schwarz -adds true -adds_iter [1] Number of iterations
*/

namespace MathLib
{
	
typedef struct {
    int ls_method;
    int ls_precond;
    std::string ls_extra_arg;
    long ls_max_iterations;
    double ls_error_tolerance;
} LIS_option;

void solveWithLis(CRS *A, double *x, double *b, LIS_option &option)
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

    sprintf(solver_options, "-i %d -p %d %s", option.ls_method, option.ls_precond, option.ls_extra_arg.c_str()); 
    sprintf(tol_option, "-tol %e -maxiter %d", option.ls_error_tolerance, option.ls_max_iterations);

    ierr = lis_matrix_set_crs(A->nonzero, A->row_ptr, A->col_idx, A->data, AA);
    ierr = lis_matrix_assemble(AA);

    // Assemble the vector, b, x
    ierr = lis_vector_duplicate(AA, &bb);
    ierr = lis_vector_duplicate(AA, &xx);
//#pragma omp parallel for
    for(size_t i=0; i < A->dimension; ++i)
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
//#pragma omp parallel for
    for(size_t i=0; i<A->dimension; ++i)
    {
        lis_vector_get_value(xx,i,&(x[i]));
    }

    // Clear memory
    lis_matrix_destroy(AA);
    lis_vector_destroy(bb);
    lis_vector_destroy(xx);
    lis_solver_destroy(solver);
    std::cout << "------------------------------------------------------------------" << std::endl;
}

}
