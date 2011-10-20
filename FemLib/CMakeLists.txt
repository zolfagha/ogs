# Source files

set( SOURCES
	Fem.h
)

include_directories(
	.
	../Base
	../GeoLib
	../MathLib
	../MeshLib
)


IF (PARALLEL_USE_MPI)
	INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF (PARALLEL_USE_MPI)

# Create the library
add_library( FEM STATIC ${SOURCES} )

target_link_libraries (
	FEM
	BASE
        MathLib
	GEO
	MeshLib
)
