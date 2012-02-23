
#include <gtest/gtest.h>

#include "MathLib/LinAlg/LinearEquations/DenseLinearEquations.h"
#include "MathLib/LinAlg/LinearEquations/SparseLinearEquations.h"
#ifdef USE_LIS
#include "MathLib/LinAlg/LinearEquations/LisInterface.h"
#endif
#ifdef USE_PETSC
#include "MathLib/LinAlg/LinearEquations/PETScLinearSolver.h"
#endif

#include "TestUtil.h"

using namespace MathLib;


struct Example1
{
    std::vector<double> mat;
    std::vector<size_t> list_dirichlet_bc_id;
    std::vector<double> list_dirichlet_bc_value;
    static const size_t dim_eqs = 9;
    std::vector<double> exH;
    RowMajorSparsity sparse;

    Example1()
    {
        double d_mat[] = {
            6.66667e-012, -1.66667e-012, 0, -1.66667e-012, -3.33333e-012, 0, 0, 0, 0, 
            -1.66667e-012, 1.33333e-011, -1.66667e-012, -3.33333e-012, -3.33333e-012, -3.33333e-012, 0, 0, 0, 
            0, -1.66667e-012, 6.66667e-012, 0, -3.33333e-012, -1.66667e-012, 0, 0, 0, 
            -1.66667e-012, -3.33333e-012, 0, 1.33333e-011, -3.33333e-012, 0, -1.66667e-012, -3.33333e-012, 0, 
            -3.33333e-012, -3.33333e-012, -3.33333e-012, -3.33333e-012, 2.66667e-011, -3.33333e-012, -3.33333e-012, -3.33333e-012, -3.33333e-012, 
            0, -3.33333e-012, -1.66667e-012, 0, -3.33333e-012, 1.33333e-011, 0, -3.33333e-012, -1.66667e-012, 
            0, 0, 0, -1.66667e-012, -3.33333e-012, 0, 6.66667e-012, -1.66667e-012, 0, 
            0, 0, 0, -3.33333e-012, -3.33333e-012, -3.33333e-012, -1.66667e-012, 1.33333e-011, -1.66667e-012, 
            0, 0, 0, 0, -3.33333e-012, -1.66667e-012, 0, -1.66667e-012, 6.66667e-012
        };
        mat.assign(d_mat, d_mat+dim_eqs*dim_eqs);
        size_t int_dirichlet_bc_id[] = {2,5,8,0,3,6};
        list_dirichlet_bc_id.assign(int_dirichlet_bc_id, int_dirichlet_bc_id+6);
        list_dirichlet_bc_value.resize(6);
        fill(list_dirichlet_bc_value.begin(), list_dirichlet_bc_value.begin()+3, .0);
        fill(list_dirichlet_bc_value.begin()+3, list_dirichlet_bc_value.end(), 1.0);
        exH.resize(9);
        for (size_t i=0; i<9; i++) {
            if (i%3==0) exH[i] = 1.0;
            if (i%3==1) exH[i] = 0.5;
            if (i%3==2) exH[i] = 0.;
        }
        sparse.resize(dim_eqs);
        for (size_t i=0; i<dim_eqs; i++) {
            for (size_t j=0; j<dim_eqs; j++) {
                if (mat[i*dim_eqs+j]!=.0)
                    sparse[i].insert(j);
            }
        }

    }
};

TEST(Math, LinearSolverDirect)
{
    // set problem
    Example1 ex1;

    // construct discrete eqs
    DenseLinearEquations eqs;
    eqs.create(ex1.sparse.size(), &ex1.sparse);

    //
    for (size_t i=0; i<ex1.dim_eqs; i++) {
        for (size_t j=0; j<ex1.dim_eqs; j++) {
            double v = ex1.mat[i*ex1.dim_eqs+j];
            eqs.addA(i, j, v);
        }
    }
    eqs.setKnownX(ex1.list_dirichlet_bc_id, ex1.list_dirichlet_bc_value);

    eqs.solve();

    ASSERT_DOUBLE_ARRAY_EQ(&ex1.exH[0], eqs.getX(), ex1.dim_eqs, 1.e-5);
}

TEST(Math, LinearSolverSparse)
{
    // set problem
    Example1 ex1;

    // construct discrete eqs
    SparseLinearEquations eqs;
    eqs.create(ex1.sparse.size(), &ex1.sparse);
    eqs.getOption().solver_type = SparseLinearEquations::SolverCG;
    eqs.getOption().precon_type = SparseLinearEquations::NONE;

    //
    for (size_t i=0; i<ex1.dim_eqs; i++) {
        for (size_t j=0; j<ex1.dim_eqs; j++) {
            double v = ex1.mat[i*ex1.dim_eqs+j];
            if (v!=.0)
                eqs.addA(i, j, v);
        }
    }
    eqs.setKnownX(ex1.list_dirichlet_bc_id, ex1.list_dirichlet_bc_value);

    eqs.solve();

    ASSERT_DOUBLE_ARRAY_EQ(&ex1.exH[0], eqs.getX(), ex1.dim_eqs, 1.e-5);
}

#ifdef USE_LIS
TEST(Math, LinearSolverLis1)
{
    // set problem
    Example1 ex1;

    // construct discrete eqs
    CRSLisSolver eqs;
    eqs.create(ex1.sparse.size(), &ex1.sparse);
    eqs.getOption().ls_method = LIS_option::CG;
    eqs.getOption().ls_precond = LIS_option::NONE;

    //
    for (size_t i=0; i<ex1.dim_eqs; i++) {
        for (size_t j=0; j<ex1.dim_eqs; j++) {
            double v = ex1.mat[i*ex1.dim_eqs+j];
            if (v!=.0)
                eqs.addA(i, j, v);
        }
    }

    eqs.setKnownX(ex1.list_dirichlet_bc_id, ex1.list_dirichlet_bc_value);

    eqs.solve();

    ASSERT_DOUBLE_ARRAY_EQ(&ex1.exH[0], eqs.getX(), ex1.dim_eqs, 1.e-5);
}

#endif

#ifdef USE_PETSC
TEST(Math, LinearSolverPETsc)
{
    // set problem
    Example1 ex1;
    std::vector<double> eqsX(ex1.dim_eqs, .0);
    {

        //
        char help[] = "Using PETSc package\n";
        int argc = 1;
        char exename[] = "test\n";
        char *argv[] = {"test\n"};
        char **tmp = argv;
        PetscInitialize(&argc, &tmp,(char *)0,help);

        int rank_p, size_p;
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank_p);
        MPI_Comm_size(PETSC_COMM_WORLD, &size_p);


        PETScLinearSolver eqs(ex1.dim_eqs);
        eqs.Init();
        eqs.set_rank_size(rank_p, size_p);

        eqs.Config(1e-10, KSPCG,PCILU);
        //eqs.Config(1e-10, KSPCG,PCILU);
        eqs.Initialize(); 

        for (int i=0; i<ex1.dim_eqs; i++)
            for (int j=0; j<ex1.dim_eqs; j++)
                eqs.addMatrixEntry(i, j, ex1.mat[i*ex1.dim_eqs+j]);

        eqs.AssembleMatrixPETSc();
        eqs.AssembleRHS_PETSc();
        eqs.AssembleUnkowns_PETSc();

        size_t nrows = ex1.list_dirichlet_bc_value.size();   
        PetscInt *rows_toberemoved = new PetscInt[nrows];

        for (size_t i=0; i<nrows; i++) 
        {
            int id = ex1.list_dirichlet_bc_id[i];
            double val = ex1.list_dirichlet_bc_value[i];
            rows_toberemoved[i] = id;
            eqs.set_bVectorEntry(id, val);         
            eqs.set_xVectorEntry(id, val);   
        }

        eqs.zeroRows_in_Matrix(nrows, rows_toberemoved);

        eqs.AssembleRHS_PETSc();
        eqs.AssembleUnkowns_PETSc();
        eqs.AssembleMatrixPETSc();

        eqs.Solver();

        eqs.UpdateSolutions(&eqsX[0], &eqsX[0]);

        delete [] rows_toberemoved;
        rows_toberemoved = NULL;  
    }

    PetscFinalize();

    ASSERT_DOUBLE_ARRAY_EQ(&ex1.exH[0], &eqsX[0], ex1.dim_eqs, 1.e-5);

}
#endif
