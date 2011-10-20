# Source files

set( SOURCES
	Mesh.h
	MeshIO.h
	MeshGenerator.h
#	Core.h
#	Edge.h
#	Element.h
#	Face.h
#	MeshEnums.h
#	Node.h
)

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
ADD_LIBRARY( MeshLib STATIC ${SOURCES} )

target_link_libraries (
	MeshLib
	Base
	GEO
        MathLib
)

SET_TARGET_PROPERTIES(MathLib PROPERTIES LINKER_LANGUAGE CXX)
