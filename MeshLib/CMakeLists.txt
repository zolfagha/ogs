# Source files
FILE(GLOB MeshLib_Core_Files Core/*.h Core/*.cpp)
SOURCE_GROUP(MeshLib\\Core FILES ${MeshLib_Core_Files})

FILE(GLOB MeshLib_IO_Files IO/*.h IO/*.cpp)
SOURCE_GROUP(MeshLib\\IO FILES ${MeshLib_IO_Files})

FILE(GLOB MeshLib_Topology_Files Topology/*.h Topology/*.cpp)
SOURCE_GROUP(MeshLib\\Topology FILES ${MeshLib_Topology_Files})

FILE(GLOB MeshLib_Tools_Files Tools/*.h Tools/*.cpp)
SOURCE_GROUP(MeshLib\\Tools FILES ${MeshLib_Tools_Files})

SET ( MeshLib_Test_Files
	Test/MeshLibTest.cpp
	)
SOURCE_GROUP(MeshLib\\Test FILES ${MeshLib_Test_Files})

ADD_LIBRARY(MeshLib STATIC 
	${MeshLib_Core_Files} 
	${MeshLib_IO_Files} 
	${MeshLib_Topology_Files} 
	${MeshLib_Tools_Files} 
	${MeshLib_Test_Files}
)

include_directories(
	.
	..
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
