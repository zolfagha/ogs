#Source files grouped by a directory
FILE(GLOB SolutionLib_Files *.h *.cpp)
SOURCE_GROUP( SolutionLib FILES ${SolutionLib_Files} )
FILE(GLOB SolutionLib_Problem_Files FemProblem/*.h FemProblem/*.cpp)
SOURCE_GROUP( SolutionLib\\FemProblem FILES ${SolutionLib_Problem_Files} )
FILE(GLOB SolutionLib_Solution_Files Solution/*.h Solution/*.cpp)
SOURCE_GROUP( SolutionLib\\Solution FILES ${SolutionLib_Solution_Files} )
FILE(GLOB SolutionLib_Tools_Files Tools/*.h Tools/*.cpp)
SOURCE_GROUP( SolutionLib\\Tools FILES ${SolutionLib_Tools_Files} )

ADD_LIBRARY( SolutionLib STATIC 
	${SolutionLib_Files}
    ${SolutionLib_Problem_Files}
    ${SolutionLib_Solution_Files}
    ${SolutionLib_Tools_Files}
)

INCLUDE_DIRECTORIES(
	.
	..
)


IF (PARALLEL_USE_MPI)
	INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF (PARALLEL_USE_MPI)

target_link_libraries (	SolutionLib
	Base
        MathLib
	GeoLib
        MeshLib
        FemLib
	NumLib
)
SET_TARGET_PROPERTIES(SolutionLib PROPERTIES LINKER_LANGUAGE CXX)
