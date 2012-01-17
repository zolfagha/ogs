#Source files grouped by a directory
SET ( ModelLib_Files
	GROUNDWATER_FLOW.h
	Test.cpp
	)
SOURCE_GROUP( ModelLib FILES ${ModelLib_Files})

ADD_LIBRARY(ModelLib STATIC ${ModelLib_Files})

INCLUDE_DIRECTORIES(
	.
	..
)


IF (PARALLEL_USE_MPI)
	INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF (PARALLEL_USE_MPI)

target_link_libraries (
	ModelLib
	BASE
        MathLib
	GeoLib
	MeshLib
	FemLib
)
