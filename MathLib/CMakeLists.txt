#Source files grouped by a directory
#./
FILE(GLOB MathLib_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
	*.h 
	*.cpp)
SOURCE_GROUP( MathLib FILES ${MathLib_Files} )
#LinAlg/
FILE(GLOB MathLib_LinAlg_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
	LinAlg/*.h 
	LinAlg/*.cpp)
SOURCE_GROUP( MathLib\\LinAlg FILES ${MathLib_LinAlg_Files} )
#LinAlg/Dense/
FILE(GLOB MathLib_LinAlg_Dense_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
	LinAlg/Dense/*.h 
	LinAlg/Dense/*.cpp)
SOURCE_GROUP( MathLib\\LinAlg\\Dense FILES ${MathLib_LinAlg_Dense_Files} )
#LinAlg/Sparse/
FILE(GLOB MathLib_LinAlg_Sparse_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
	LinAlg/Sparse/*.h 
	LinAlg/Sparse/*.cpp)
IF (NOT USE_EIGEN)
	LIST(REMOVE_ITEM MathLib_LinAlg_Sparse_Files 
		LinAlg/Sparse/EigenInterface.h)
ENDIF ()
SOURCE_GROUP( MathLib\\LinAlg\\Sparse FILES ${MathLib_LinAlg_Sparse_Files} )
#LinAlg/Solvers/
FILE(GLOB MathLib_LinAlg_Solvers_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} LinAlg/Solvers/*.h LinAlg/Solvers/*.cpp)
IF (NOT USE_BLAS_LAPACK)
	LIST(REMOVE_ITEM MathLib_LinAlg_Solvers_Files 
	        LinAlg/Solvers/BiCGStab.h
	        LinAlg/Solvers/BiCGStab.cpp
	        LinAlg/Solvers/CG.h
	        LinAlg/Solvers/CG.cpp
	)
ENDIF ()
SOURCE_GROUP( MathLib\\LinAlg\\Solvers FILES ${MathLib_LinAlg_Solvers_Files} )
#LinAlg/Preconditioner/
FILE(GLOB MathLib_LinAlg_Preconditioner_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
	LinAlg/Preconditioner/*.h 
	LinAlg/Preconditioner/*.cpp)
SOURCE_GROUP( MathLib\\LinAlg\\Preconditioner FILES ${MathLib_LinAlg_Preconditioner_Files} )
#LinAlg/LinearEquations/
FILE(GLOB MathLib_LinAlg_LinearEquations_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
	LinAlg/LinearEquations/*.h 
	LinAlg/LinearEquations/*.cpp)
IF (NOT USE_BLAS_LAPACK)
	LIST(REMOVE_ITEM MathLib_LinAlg_LinearEquations_Files 
		LinAlg/LinearEquations/SparseLinearEquations.h
		LinAlg/LinearEquations/SparseLinearEquations.cpp
	)
ENDIF ()
IF (NOT USE_LIS)
	LIST(REMOVE_ITEM MathLib_LinAlg_LinearEquations_Files 
		LinAlg/LinearEquations/LisInterface.h
		LinAlg/LinearEquations/LisInterface.cpp
	)
ENDIF ()
IF(NOT (USE_MPI AND USE_LIS))
	LIST(REMOVE_ITEM MathLib_LinAlg_LinearEquations_Files 
		LinAlg/LinearEquations/LisMPILinearEquation.h
		LinAlg/LinearEquations/LisMPILinearEquation.cpp
	)
ENDIF()

IF (NOT USE_PETSC)
	LIST(REMOVE_ITEM MathLib_LinAlg_LinearEquations_Files 
		LinAlg/LinearEquations/PETScLinearSolver.h
		LinAlg/LinearEquations/PETScLinearSolver.cpp
	)
ENDIF ()
SOURCE_GROUP( MathLib\\LinAlg\\LinearEquations FILES ${MathLib_LinAlg_LinearEquations_Files})
#Nonlinear/
FILE(GLOB MathLib_Nonlinear_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
	Nonlinear/*.h 
	Nonlinear/*.cpp)
SOURCE_GROUP( MathLib\\Nonlinear FILES ${MathLib_Nonlinear_Files})
#Integration/
FILE(GLOB MathLib_Integration_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
	Integration/*.h 
	Integration/*.cpp)
SOURCE_GROUP( MathLib\\Integration FILES ${MathLib_Integration_Files})
#ODE/
FILE(GLOB MathLib_ODE_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
	ODE/*.h 
	ODE/*.cpp)
SOURCE_GROUP( MathLib\\ODE FILES ${MathLib_ODE_Files})
#Function/
FILE(GLOB MathLib_Function_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
	Function/*.h 
	Function/*.cpp)
SOURCE_GROUP( MathLib\\Function FILES ${MathLib_Function_Files})
#Parameter/
FILE(GLOB MathLib_Parameter_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
    Parameter/*.h 
    Parameter/*.cpp)
SOURCE_GROUP( MathLib\\Parameter FILES ${MathLib_Parameter_Files})
#Coupling/
FILE(GLOB MathLib_Coupling_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
    Coupling/*.h 
    Coupling/*.cpp)
SOURCE_GROUP( MathLib\\Coupling FILES ${MathLib_Coupling_Files})
#Coupling/Algorithm
FILE(GLOB MathLib_Coupling_Algorithm_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
    Coupling/Algorithm/*.h 
    Coupling/Algorithm/*.cpp)
SOURCE_GROUP( MathLib\\Coupling\\Algorithm FILES ${MathLib_Coupling_Algorithm_Files})
#SystemOfEquations/
FILE(GLOB MathLib_SystemOfEquations_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
    SystemOfEquations/*.h 
    SystemOfEquations/*.cpp)
SOURCE_GROUP( MathLib\\SystemOfEquations FILES ${MathLib_SystemOfEquations_Files})

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
    ${MathLib_Parameter_Files}
    ${MathLib_Coupling_Files}
    ${MathLib_Coupling_Algorithm_Files}
    ${MathLib_SystemOfEquations_Files}
)

#Include
INCLUDE_DIRECTORIES (
	.
	..
)
IF (USE_MPI)
	INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF (USE_MPI)

IF(USE_LIS)
	INCLUDE_DIRECTORIES( ${LIS_INCLUDE_DIR} )
ENDIF(USE_LIS)

SET_TARGET_PROPERTIES(MathLib PROPERTIES LINKER_LANGUAGE CXX)

