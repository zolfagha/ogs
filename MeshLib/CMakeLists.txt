# Source files
SET ( MeshLib_Files
	Mesh.h
	MeshIO.h
	MeshGenerator.h
	)
SOURCE_GROUP(MeshLib FILES ${MeshLib_Files})

ADD_LIBRARY(MeshLib STATIC ${MeshLib_Files})

include_directories(
	.
	../Base
	../GeoLib
	../MathLib
)


IF (PARALLEL_USE_MPI)
	INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF (PARALLEL_USE_MPI)

# Create the library

target_link_libraries (
	MeshLib
	Base
	GeoLib
        MathLib
)

SET_TARGET_PROPERTIES(MathLib PROPERTIES LINKER_LANGUAGE CXX)
