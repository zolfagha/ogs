#Source files grouped by a directory
FILE(GLOB SolutionLib_Files *.h *.cpp)
SOURCE_GROUP( SolutionLib FILES ${SolutionLib_Files} )

ADD_LIBRARY( SolutionLib STATIC 
	${SolutionLib_Files}
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
SET_TARGET_PROPERTIES(Base PROPERTIES LINKER_LANGUAGE CXX)
