# Source files
FILE(GLOB MeshLib_Core_Files MeshLib/Core/*.h MeshLib/Core/*.cpp)
SOURCE_GROUP(ogsMesh\\Core FILES ${MeshLib_Core_Files})

FILE(GLOB MeshLib_IO_Files MeshLib/IO/*.h MeshLib/IO/*.cpp)
SOURCE_GROUP(ogsMesh\\IO FILES ${MeshLib_IO_Files})

FILE(GLOB MeshLib_Topology_Files MeshLib/Topology/*.h MeshLib/Topology/*.cpp)
SOURCE_GROUP(ogsMesh\\Topology FILES ${MeshLib_Topology_Files})

FILE(GLOB MeshLib_Tools_Files MeshLib/Tools/*.h MeshLib/Tools/*.cpp)
SOURCE_GROUP(ogsMesh\\Tools FILES ${MeshLib_Tools_Files})

#SET ( MeshLib_Test_Files
#	Test/MeshLibTest.cpp
#	)
#SOURCE_GROUP(ogsNumerics\\MeshLib\\Test FILES ${MeshLib_Test_Files})

#
LIST(APPEND MeshLib_Files
	${MeshLib_Core_Files} 
	${MeshLib_IO_Files} 
	${MeshLib_Topology_Files} 
	${MeshLib_Tools_Files} 
	${MeshLib_Test_Files}
)

#ADD_LIBRARY(ogsMesh STATIC 
#	${MeshLib_Core_Files} 
#	${MeshLib_IO_Files} 
#	${MeshLib_Topology_Files} 
#	${MeshLib_Tools_Files} 
#	${MeshLib_Test_Files}
#)
#
#INCLUDE_DIRECTORIES(
#	.
#	..
#	../ogsCore
#)
#

#IF (PARALLEL_USE_MPI)
#	INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
#ENDIF (PARALLEL_USE_MPI)
#
## Create the library
#
#target_link_libraries (
#	MeshLib
#	BaseLib
#	GeoLib
#    MathLib
#)
#
#SET_TARGET_PROPERTIES(ogsMesh PROPERTIES LINKER_LANGUAGE CXX)

