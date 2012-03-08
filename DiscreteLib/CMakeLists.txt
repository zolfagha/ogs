#Source files grouped by a directory
FILE(GLOB DiscreteLib_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h *.cpp)
IF(NOT (USE_LIS AND USE_MPI))
	LIST(REMOVE_ITEM DiscreteLib_Files LisMPIDiscreteVector.h LisMPIDiscreteVector.cpp)
ENDIF()
SOURCE_GROUP( DiscreteLib FILES ${DiscreteLib_Files} )

FILE(GLOB DiscreteLib_OGS5_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ogs5/*.h ogs5/*.cpp)
SOURCE_GROUP( DiscreteLib\\OGS5 FILES ${DiscreteLib_OGS5_Files} )


ADD_LIBRARY( DiscreteLib STATIC 
	${DiscreteLib_Files} 
	${DiscreteLib_OGS5_Files} 
)

INCLUDE_DIRECTORIES(
	.
	..
)

IF (USE_MPI)
	INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF (USE_MPI)

target_link_libraries (
	DiscreteLib
	Base
        MathLib
        MeshLib
)

SET_TARGET_PROPERTIES(DiscreteLib PROPERTIES LINKER_LANGUAGE CXX)
