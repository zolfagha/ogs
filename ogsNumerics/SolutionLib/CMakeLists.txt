#Source files grouped by a directory
FILE(GLOB SolutionLib_Files SolutionLib/*.h SolutionLib/*.cpp)
SOURCE_GROUP( ogsNumerics\\SolutionLib FILES ${SolutionLib_Files} )
FILE(GLOB SolutionLib_Problem_Files SolutionLib/FemProblem/*.h SolutionLib/FemProblem/*.cpp)
SOURCE_GROUP( ogsNumerics\\SolutionLib\\FemProblem FILES ${SolutionLib_Problem_Files} )
FILE(GLOB SolutionLib_Solution_Files SolutionLib/Solution/*.h SolutionLib/Solution/*.cpp)
SOURCE_GROUP( ogsNumerics\\SolutionLib\\Solution FILES ${SolutionLib_Solution_Files} )
FILE(GLOB SolutionLib_Tools_Files SolutionLib/Tools/*.h SolutionLib/Tools/*.cpp)
SOURCE_GROUP( ogsNumerics\\SolutionLib\\Tools FILES ${SolutionLib_Tools_Files} )

#
LIST(APPEND SolutionLib_Files
	${SolutionLib_Files}
    ${SolutionLib_Problem_Files}
    ${SolutionLib_Solution_Files}
    ${SolutionLib_Tools_Files}
)

#ADD_LIBRARY( SolutionLib STATIC 
#	${SolutionLib_Files}
#    ${SolutionLib_Problem_Files}
#    ${SolutionLib_Solution_Files}
#    ${SolutionLib_Tools_Files}
#)
#
#INCLUDE_DIRECTORIES(
#	.
#	..
#	../../ogsCore
#)
#
#
#IF (PARALLEL_USE_MPI)
#	INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
#ENDIF (PARALLEL_USE_MPI)
#
#target_link_libraries (	SolutionLib
#	BaseLib
#        MathLib
#	GeoLib
#        MeshLib
#        FemLib
#	NumLib
#)
#SET_TARGET_PROPERTIES(SolutionLib PROPERTIES LINKER_LANGUAGE CXX)
#