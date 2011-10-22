#Source files grouped by a directory
SET ( FemLib_Files
	Fem.h
	)
SOURCE_GROUP( FemLib FILES ${FemLib_Files})

ADD_LIBRARY(FemLib STATIC ${FemLib_Files})

INCLUDE_DIRECTORIES(
	.
	../Base
	../GeoLib
	../MathLib
	../MeshLib
)


IF (PARALLEL_USE_MPI)
	INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF (PARALLEL_USE_MPI)

target_link_libraries (
	FemLib
	BASE
        MathLib
	GEO
	MeshLib
)
