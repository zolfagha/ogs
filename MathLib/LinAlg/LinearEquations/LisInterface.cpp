
#include "LisInterface.h"

#include <iostream>
#ifdef _OPENMP
#include <omp.h>
#endif


namespace MathLib
{

void CRSLisSolver::initialize()
{
    int argc=0;
    char **argv;
    lis_initialize(&argc, &argv);
}

void CRSLisSolver::finalize()
{
    lis_finalize();
}

void CRSLisSolver::setOption(const Base::Options &option)
{
    throw std::exception("LisSolver::setOption() is not implemented.");
}

void CRSLisSolver::solveEqs(CRSMatrix<double, signed> *A, double *b, double *x)
{
    long dimension = static_cast<long>(A->getNRows());

    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "*** LIS solver computation" << std::endl;

    // Creating a matrix.
    LIS_MATRIX AA;
    LIS_VECTOR bb,xx;
    LIS_SOLVER solver;
#ifndef USE_MPI
    int ierr = lis_matrix_create(0, &AA);
    ierr = lis_matrix_set_type(AA, LIS_MATRIX_CRS);
    ierr = lis_matrix_set_size(AA, 0, dimension);
#else
    lis_matrix_create(MPI_COMM_WORLD, &AA);
    lis_matrix_set_size(AA,dimension,0);
//    lis_matrix_get_size(AA, &_local_dim, &_global_dim);
#endif

    // Matrix solver and Precondition can be handled better way.
    const size_t MAX_ZEILE = 512;
    char solver_options[MAX_ZEILE], tol_option[MAX_ZEILE];

    int nthreads = omp_get_num_threads();
    //omp_set_num_threads (nthreads);

    sprintf(solver_options, "-i %d -p %d %s", _option.ls_method, _option.ls_precond, _option.ls_extra_arg.c_str()); 
    sprintf(tol_option, "-tol %e -maxiter %d -omp_num_threads %d -initx_zeros 0", _option.ls_error_tolerance, _option.ls_max_iterations, nthreads);

    ierr = lis_matrix_set_crs(A->getNNZ(), (int*)A->getRowPtrArray(), (int*)A->getColIdxArray(), (double*)A->getEntryArray(), AA);
    ierr = lis_matrix_assemble(AA);

    // Assemble the vector, b, x
    ierr = lis_vector_duplicate(AA, &bb);
    ierr = lis_vector_duplicate(AA, &xx);
    ierr = lis_vector_scatter(x, xx);
    ierr = lis_vector_scatter(b, bb);
    //#pragma omp parallel for
    //for (long i=0; i < dimension; ++i)
    //{
    //    ierr = lis_vector_set_value(LIS_INS_VALUE, i, x[i], xx);
    //    ierr = lis_vector_set_value(LIS_INS_VALUE, i, b[i], bb);
    //}

    // Create solver
    ierr = lis_solver_create(&solver);

    ierr = lis_solver_set_option(solver_options, solver);
    ierr = lis_solver_set_option(tol_option, solver);
    ierr = lis_solver_set_option("-print mem", solver);
    
    ierr = lis_solve(AA, bb, xx, solver);
    int iter = 0;
    double resid = 0.0;
    ierr = lis_solver_get_iters(solver, &iter);
    ierr = lis_solver_get_residualnorm(solver, &resid);
    printf("\t iteration: %d/%d\n", iter, _option.ls_max_iterations);
    printf("\t residuals: %e\n", resid);
    //	lis_vector_print(xx);
    //	lis_vector_print(bb);

    // Update the solution (answer) into the x vector
    ierr = lis_vector_gather(xx, x);

    //#pragma omp parallel for
    //for(long i=0; i<dimension; ++i)
    //{
    //    lis_vector_get_value(xx,i,&(x[i]));
    //}

    // Clear memory
    //lis_matrix_destroy(AA);
    lis_vector_destroy(bb);
    lis_vector_destroy(xx);
    lis_solver_destroy(solver);
    std::cout << "------------------------------------------------------------------" << std::endl;
}

#if 0
void LisSolver::initialize()
{
    int argc=0;
    char **argv;
    lis_initialize(&argc, &argv);
}

void LisSolver::finalize()
{
    lis_finalize();
}

LisSolver::~LisSolver()
{
    lis_matrix_destroy(_A);
    lis_vector_destroy(_b);
    lis_vector_destroy(_x);
}

void LisSolver::create(size_t length, RowMajorSparsity *sparsity)
{
    int n = static_cast<int>(length);
    int is = 0;
    int ie = n;
#ifndef USE_MPI
    lis_matrix_create(0, &_A);
    lis_matrix_set_size(_A,0,n);
    lis_vector_create(0, &_b);
    lis_vector_create(0, &_x);
    lis_vector_set_size(_b, 0, n);
    lis_vector_set_size(_x, 0, n);
#else
    lis_matrix_create(MPI_COMM_WORLD, &_A);
    lis_matrix_set_size(_A,0,n);
    lis_matrix_get_size(_A, &_local_dim, &_global_dim);
    lis_vector_create(MPI_COMM_WORLD, &_b);
    lis_vector_create(MPI_COMM_WORLD, &_x);
    lis_vector_get_range(_b, &is, &ie);
#endif

    reset();
}

void LisSolver::setOption(const Base::Options &option)
{
    throw std::exception("LisSolver::setOption() is not implemented.");
}

void LisSolver::reset()
{
    lis_vector_set_all(0., _b);
    lis_vector_set_all(0., _x);
}

double LisSolver::getA(size_t rowId, size_t colId)
{
    throw std::exception("not implemented.");
}

void LisSolver::setA(size_t rowId, size_t colId, double v)
{
    lis_matrix_set_value(LIS_INS_VALUE, rowId, colId, v, _A);
}

void LisSolver::addA(size_t rowId, size_t colId, double v)
{
    lis_matrix_set_value(LIS_ADD_VALUE, rowId, colId, v, _A);
}

void LisSolver::addA(std::vector<size_t> &vec_row_pos, std::vector<size_t> &vec_col_pos, MathLib::Matrix<double> &sub_matrix, double fkt)
{
    for (size_t i=0; i<vec_row_pos.size(); i++) {
        const size_t rowId = vec_row_pos[i];
        for (size_t j=0; j<vec_col_pos.size(); j++) {
            const size_t colId = vec_col_pos[j];
            addA(rowId, colId, fkt*sub_matrix(i,j));
        }
    }
}

void LisSolver::addA(std::vector<size_t> &vec_pos, MathLib::Matrix<double> &sub_matrix, double fkt)
{
    addA(vec_pos, vec_pos, sub_matrix, fkt);
}

double LisSolver::getRHS(size_t rowId)
{
    double v;
    lis_vector_get_value(_b, rowId, &v);
    return v;
}

double* LisSolver::getRHS()
{
    return &_tmp_b[0];
}

void LisSolver::setRHS(size_t rowId, double v)
{
    lis_vector_set_value(LIS_INS_VALUE, rowId, v, _b);
}

void LisSolver::addRHS(size_t rowId, double v)
{
    lis_vector_set_value(LIS_ADD_VALUE, rowId, v, _b);
}

void LisSolver::addRHS(std::vector<size_t> &vec_pos, double *sub_vector, double fkt)
{
    for (size_t i=0; i<vec_pos.size(); i++) {
        const size_t rowId = vec_pos[i];
        addRHS(rowId, sub_vector[i]*fkt);
    }
}

double* LisSolver::getX()
{
    return &_tmp_x[0];
}

void LisSolver::setKnownX(size_t row_id, double x)
{

}

void LisSolver::setKnownX(const std::vector<size_t> &vec_id, const std::vector<double> &vec_x)
{

}

void LisSolver::solve()
{
    int iter = 0;
    double resid = 0.0;

    //assemble matrix
    lis_matrix_set_type(_A, LIS_MATRIX_CRS);
    lis_matrix_assemble(_A);




    // Matrix solver and Precondition can be handled better way.
    const size_t MAX_ZEILE = 512;
    char solver_options[MAX_ZEILE], tol_option[MAX_ZEILE];

    int nthreads = omp_get_num_threads();
    //omp_set_num_threads (nthreads);

    sprintf(solver_options, "-i %d -p %d %s", _option.ls_method, _option.ls_precond, _option.ls_extra_arg.c_str()); 
    sprintf(tol_option, "-tol %e -maxiter %d -omp_num_threads %d -initx_zeros 0", _option.ls_error_tolerance, _option.ls_max_iterations, nthreads);



    // Create solver
    LIS_SOLVER solver;
    lis_solver_create(&solver);
    lis_solver_set_option(solver_options, solver);
    lis_solver_set_option(tol_option, solver);
    lis_solver_set_option("-print mem", solver);

    // solve
    lis_solve(_A, _b, _x, solver);
    _tmp_x.resize(_global_dim);
    lis_vector_gather(_x, &_tmp_x[0]);

    //
    lis_solver_get_iters(solver, &iter);
    lis_solver_get_residualnorm(solver, &resid);
    printf("\t iteration: %d/%d\n", iter, _option.ls_max_iterations);
    printf("\t residuals: %e\n", resid);

    lis_solver_destroy(solver);
}
#endif

}
