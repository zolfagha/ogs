
#include <gtest/gtest.h>

#include "MathLib/LinAlg/Dense/Matrix.h"
#include "MathLib/LinAlg/LinearEquations/LinearEquationsFactory.h"

#include "MathLib/LinAlg/LinearEquations/SparseLinearEquations.h"
#ifdef USE_PETSC
#include "MathLib/LinAlg/LinearEquations/PETScLinearSolver.h"
#endif

#include "TestUtil.h"

using namespace MathLib;

TEST(Math, Matrix_transposeAndMultiply)
{
    MathLib::Matrix<double> matA(1,2);
    MathLib::Matrix<double> matB(1,2);
    MathLib::Matrix<double> matC(2,2);

    matA(0,0) = 1.0;
    matA(0,1) = 2.0;
    matB(0,0) = 3.0;
    matB(0,1) = 4.0;
    matC = .0;
    matA.transposeAndMultiply(matB, matC);

    ASSERT_EQ(matC(0,0), 3.0);
    ASSERT_EQ(matC(0,1), 4.0);
    ASSERT_EQ(matC(1,0), 6.0);
    ASSERT_EQ(matC(1,1), 8.0);
}

TEST(Math, SparseLinearEQS)
{
    // set problem
    const size_t dim_eqs = 9;
    size_t int_dirichlet_bc_id[] = {2,5,8,0,3,6};
    std::vector<size_t> list_dirichlet_bc_id(int_dirichlet_bc_id, int_dirichlet_bc_id+6);
    std::vector<double> list_dirichlet_bc_value(6);
    fill(list_dirichlet_bc_value.begin(), list_dirichlet_bc_value.begin()+3, .0);
    fill(list_dirichlet_bc_value.begin()+3, list_dirichlet_bc_value.end(), 1.0);

    double mat[] = {
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

    //
    RowMajorSparsity sparse;
    sparse.resize(dim_eqs);
    for (size_t i=0; i<dim_eqs; i++) {
        for (size_t j=0; j<dim_eqs; j++) {
            if (mat[i*dim_eqs+j]!=.0)
                sparse[i].insert(j);
        }
    }

    // construct discrete eqs
    SparseLinearEquations eqs;
    eqs.create(sparse.size(), &sparse);
    eqs.getOption().solver_type = SparseLinearEquations::SolverCG;
    eqs.getOption().precon_type = SparseLinearEquations::NONE;

    //
    for (size_t i=0; i<dim_eqs; i++) {
        for (size_t j=0; j<dim_eqs; j++) {
            double v = mat[i*dim_eqs+j];
            if (v!=.0)
                eqs.addA(i, j, v);
        }
    }

    eqs.setKnownX(list_dirichlet_bc_id, list_dirichlet_bc_value);

    eqs.solve();

    std::vector<double> exH(9);
    for (size_t i=0; i<9; i++) {
        if (i%3==0) exH[i] = 1.0;
        if (i%3==1) exH[i] = 0.5;
        if (i%3==2) exH[i] = 0.;
    }
    ASSERT_DOUBLE_ARRAY_EQ(&exH[0], eqs.getX(), dim_eqs, 1.e-5);
}

#ifdef USE_PETSC
TEST(Math, PETSC)
{
    // set problem
    const size_t dim_eqs = 9;

    std::vector<double> eqsX(dim_eqs, .0);

    size_t list_dirichlet_bc_id[] = {2,5,8,0,3,6};
    std::vector<double> list_dirichlet_bc_value(6);
    fill(list_dirichlet_bc_value.begin(), list_dirichlet_bc_value.begin()+3, .0);
    fill(list_dirichlet_bc_value.begin()+3, list_dirichlet_bc_value.end(), 1.0);

    double mat[] = {
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


    PETScLinearSolver eqs(dim_eqs);
    eqs.Init();
#if 1
    eqs.set_rank_size(rank_p, size_p);

    eqs.Config(1e-10, KSPCG,PCILU);
    //eqs.Config(1e-10, KSPCG,PCILU);
    eqs.Initialize(); 

    for (int i=0; i<dim_eqs; i++)
        for (int j=0; j<dim_eqs; j++)
            eqs.addMatrixEntry(i, j, mat[i*dim_eqs+j]);

    eqs.AssembleMatrixPETSc();
    eqs.AssembleRHS_PETSc();
    eqs.AssembleUnkowns_PETSc();

    int nrows = (int)list_dirichlet_bc_value.size();   
    PetscInt *rows_toberemoved = new PetscInt[nrows];

    for (size_t i=0; i<nrows; i++) 
    {
        int id = list_dirichlet_bc_id[i];
        double val = list_dirichlet_bc_value[i];
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
#endif
    }

    PetscFinalize();

    std::vector<double> exH(9);
    for (size_t i=0; i<9; i++) {
        if (i%3==0) exH[i] = 1.0;
        if (i%3==1) exH[i] = 0.5;
        if (i%3==2) exH[i] = 0.;
    }
    ASSERT_DOUBLE_ARRAY_EQ(&exH[0], &eqsX[0], dim_eqs, 1.e-5);

}
#endif
