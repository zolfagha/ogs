#Source files grouped by a directory
SET ( MathLib_Files
	AnalyticalGeometry.h        
	LinearInterpolation.h  
	MathTools.h
	Vector3.h
	EarClippingTriangulation.h  
	max.h
        sparse.h
	vector_io.h
	AnalyticalGeometry.cpp        
	LinearInterpolation.cpp  
	MathTools.cpp
	EarClippingTriangulation.cpp  
	)
SOURCE_GROUP( MathLib FILES ${MathLib_Files})

SET ( MathLib_LinAlg_Files
	LinAlg/MatrixBase.h
	LinAlg/VectorNorms.h
	LinAlg/Dense/Matrix.h
        LinAlg/Preconditioner/generateDiagPrecond.h
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
	LinAlg/Solvers/LisInterface.h
	LinAlg/Sparse/amuxCRS.h
        LinAlg/Sparse/CRSMatrix.h
        LinAlg/Sparse/CRSMatrixPThreads.h
        LinAlg/Sparse/CRSMatrixOpenMP.h
        LinAlg/Sparse/CRSSymMatrix.h
        LinAlg/Sparse/SparseMatrixBase.h
        LinAlg/Sparse/EigenInterface.h
        LinAlg/Sparse/amuxCRS.cpp
        LinAlg/Solvers/BiCGStab.cpp
        LinAlg/Solvers/CG.cpp
        LinAlg/Solvers/GMRes.cpp
	LinAlg/Solvers/GaussAlgorithm.cpp
        LinAlg/Solvers/TriangularSolve.cpp
	LinAlg/Preconditioner/generateDiagPrecond.cpp
	)
SOURCE_GROUP( MathLib\\LinAlg FILES ${MathLib_LinAlg_Files})

# Create the library
ADD_LIBRARY( MathLib STATIC ${MathLib_Files} ${MathLib_LinAlg_Files} )

#Include
INCLUDE_DIRECTORIES (
	.
	../Base
	../GeoLib
)

SET_TARGET_PROPERTIES(MathLib PROPERTIES LINKER_LANGUAGE CXX)

