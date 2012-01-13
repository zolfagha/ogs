# Source files
SET ( MeshLib_Files
	Mesh.h
	Node.h
	Element.h
	ElementFactory.h
	)
SOURCE_GROUP(MeshLib FILES ${MeshLib_Files})

SET ( MeshLib_IO_Files
	IO/MeshIOOGSAscii.h
	IO/MeshIOLegacyVtk.h
	)
SOURCE_GROUP(MeshLib\\IO FILES ${MeshLib_IO_Files})

SET ( MeshLib_Tools_Files
	Tools/MeshGenerator.h
	)
SOURCE_GROUP(MeshLib\\Tools FILES ${MeshLib_Tools_Files})

ADD_LIBRARY(MeshLib STATIC ${MeshLib_Files} ${MeshLib_IO_Files} ${MeshLib_Tools_Files})

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
