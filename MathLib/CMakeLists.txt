#Source files grouped by a directory
SET ( MathLib_Files
	LinearInterpolation.h  
	MathTools.h
	Vector.h
	max.h
	sparse.h
	vector_io.h
	LinearInterpolation.cpp  
	MathTools.cpp
)
SOURCE_GROUP( MathLib FILES ${MathLib_Files})

SET ( MathLib_LinAlg_Files
	LinAlg/MatrixBase.h
	LinAlg/VectorNorms.h
)
SOURCE_GROUP( MathLib\\LinAlg FILES ${MathLib_LinAlg_Files})

SET ( MathLib_LinAlg_Dense_Files
	LinAlg/Dense/Matrix.h
)
SOURCE_GROUP( MathLib\\LinAlg\\Dense FILES ${MathLib_LinAlg_Dense_Files})

SET ( MathLib_LinAlg_Sparse_Files
	LinAlg/Sparse/amuxCRS.h
        LinAlg/Sparse/CRSMatrix.h
        LinAlg/Sparse/CRSMatrixPThreads.h
        LinAlg/Sparse/CRSMatrixOpenMP.h
        LinAlg/Sparse/CRSSymMatrix.h
        LinAlg/Sparse/SparseMatrixBase.h
        LinAlg/Sparse/amuxCRS.cpp
        LinAlg/Sparse/SparseTableCRS.h
)
IF (USE_EIGEN)
SET ( MathLib_LinAlg_Sparse_Files ${MathLib_LinAlg_Sparse_Files}
	        LinAlg/Sparse/EigenInterface.h
)
ENDIF (USE_EIGEN)
SOURCE_GROUP( MathLib\\LinAlg\\Sparse FILES ${MathLib_LinAlg_Sparse_Files})

SET ( MathLib_LinAlg_Solvers_Files
        LinAlg/Solvers/LinearSolver.h
        LinAlg/Solvers/DirectLinearSolver.h
        LinAlg/Solvers/DenseDirectLinearSolver.h
        LinAlg/Solvers/GaussAlgorithm.h
        LinAlg/Solvers/TriangularSolve.h
        LinAlg/Solvers/IterativeLinearSolver.h
        LinAlg/Solvers/solver.h
        LinAlg/Solvers/BiCGStab.h
        LinAlg/Solvers/CG.h
        LinAlg/Solvers/GMRes.h
        LinAlg/Solvers/BiCGStab.cpp
        LinAlg/Solvers/CG.cpp
        LinAlg/Solvers/CGParallel.cpp
        LinAlg/Solvers/GMRes.cpp
		LinAlg/Solvers/GaussAlgorithm.cpp
        LinAlg/Solvers/TriangularSolve.cpp
)
SOURCE_GROUP( MathLib\\LinAlg\\Solvers FILES ${MathLib_LinAlg_Solvers_Files})

SET ( MathLib_LinAlg_Preconditioner_Files
		LinAlg/Preconditioner/generateDiagPrecond.h
		LinAlg/Preconditioner/generateDiagPrecond.cpp
)
SOURCE_GROUP( MathLib\\LinAlg\\Preconditioner FILES ${MathLib_LinAlg_Preconditioner_Files})

SET ( MathLib_LinAlg_LinearEquations_Files
		LinAlg/LinearEquations/ILinearEquations.h
		LinAlg/LinearEquations/DenseLinearEquations.h
		LinAlg/LinearEquations/DenseLinearEquations.cpp
		LinAlg/LinearEquations/SparseLinearEquations.h
		LinAlg/LinearEquations/SparseLinearEquations.cpp
		LinAlg/LinearEquations/PETScInterface.h
		LinAlg/LinearEquations/PARDISOInterface.h
		LinAlg/LinearEquations/LinearEquationsFactory.h
)
IF (USE_LIS)
	SET ( MathLib_LinAlg_LinearEquations_Files ${MathLib_LinAlg_LinearEquations_Files}
		LinAlg/LinearEquations/LisInterface.h
		LinAlg/LinearEquations/LisInterface.cpp
)
ENDIF (USE_LIS)
SOURCE_GROUP( MathLib\\LinAlg\\LinearEquations FILES ${MathLib_LinAlg_LinearEquations_Files})

SET ( MathLib_Nonlinear_Files
        Nonlinear/Picard.h
        Nonlinear/NewtonRaphson.h
        Nonlinear/Convergence.h
)
SOURCE_GROUP( MathLib\\Nonlinear FILES ${MathLib_Nonlinear_Files})

SET ( MathLib_Integration_Files
        Integration/GaussLegendre.h
        Integration/NewtonCotes.h
)
SOURCE_GROUP( MathLib\\Integration FILES ${MathLib_Integration_Files})

SET ( MathLib_ODE_Files
        ODE/ForwardEuler.h
        ODE/RungeKutta4.h
)
SOURCE_GROUP( MathLib\\ODE FILES ${MathLib_ODE_Files})

SET ( MathLib_Function_Files
        Function/Function.h
)
SOURCE_GROUP( MathLib\\Function FILES ${MathLib_Function_Files})


# Create the library
ADD_LIBRARY( MathLib STATIC 
	${MathLib_Files} 
	${MathLib_LinAlg_Files} 
	${MathLib_LinAlg_Dense_Files}
	${MathLib_LinAlg_Sparse_Files}
	${MathLib_LinAlg_Solvers_Files}
	${MathLib_LinAlg_Preconditioner_Files}
	${MathLib_LinAlg_LinearEquations_Files}
	${MathLib_Nonlinear_Files}
	${MathLib_Integration_Files}
	${MathLib_ODE_Files}
	${MathLib_Function_Files}
)

#Include
INCLUDE_DIRECTORIES (
	.
	..
)

IF(USE_LIS)
	INCLUDE_DIRECTORIES( ${LIS_INCLUDE_DIR} )
ENDIF(USE_LIS)

SET_TARGET_PROPERTIES(MathLib PROPERTIES LINKER_LANGUAGE CXX)

