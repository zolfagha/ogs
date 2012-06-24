#Source files grouped by a directory
#.
FILE(GLOB FemLib_Files *.h *.cpp)
SOURCE_GROUP( FemLib FILES ${FemLib_Files})
#Core
FILE(GLOB FemLib_Core_Files Core/*.h Core/*.cpp)
SOURCE_GROUP( FemLib\\Core FILES ${FemLib_Core_Files})
#Core/Element
FILE(GLOB FemLib_Core_Element_Files Core/Element/*.h Core/Element/*.cpp)
SOURCE_GROUP( FemLib\\Core\\Element FILES ${FemLib_Core_Element_Files})
#Core/CoordinatesMapping
FILE(GLOB FemLib_Core_CoordinatesMapping_Files Core/CoordinatesMapping/*.h Core/CoordinatesMapping/*.cpp)
SOURCE_GROUP( FemLib\\Core\\CoordinatesMapping FILES ${FemLib_Core_CoordinatesMapping_Files})
#Core/Shape
FILE(GLOB FemLib_Core_Shape_Files Core/Shape/*.h Core/Shape/*.cpp)
SOURCE_GROUP( FemLib\\Core\\Shape FILES ${FemLib_Core_Shape_Files})
#Core/Extrapolation
FILE(GLOB FemLib_Core_Extrapolation_Files Core/Extrapolation/*.h Core/Extrapolation/*.cpp)
SOURCE_GROUP( FemLib\\Core\\Extrapolation FILES ${FemLib_Core_Extrapolation_Files})
#Core/Integration
FILE(GLOB FemLib_Core_Integration_Files Core/Integration/*.h Core/Integration/*.cpp)
SOURCE_GROUP( FemLib\\Core\\Integration FILES ${FemLib_Core_Integration_Files})
#BC
FILE(GLOB FemLib_BC_Files BC/*.h BC/*.cpp)
SOURCE_GROUP( FemLib\\BC FILES ${FemLib_BC_Files})
#Function
FILE(GLOB FemLib_Function_Files Function/*.h Function/*.cpp)
SOURCE_GROUP( FemLib\\Function FILES ${FemLib_Function_Files})
#Post
FILE(GLOB FemLib_Post_Files Post/*.h Post/*.cpp)
SOURCE_GROUP( FemLib\\Post FILES ${FemLib_Post_Files})
#Tools
FILE(GLOB FemLib_Tools_Files Tools/*.h Tools/*.cpp)
SOURCE_GROUP( FemLib\\Tools FILES ${FemLib_Tools_Files})


ADD_LIBRARY(FemLib STATIC
	${FemLib_Files}
	${FemLib_Core_Files}
	${FemLib_Core_Element_Files}
    ${FemLib_Core_Shape_Files}
    ${FemLib_Core_Extrapolation_Files}
    ${FemLib_Core_Integration_Files}
    ${FemLib_Core_CoordinatesMapping_Files}
	${FemLib_BC_Files}
	${FemLib_Function_Files}
	${FemLib_Post_Files}
    ${FemLib_Tools_Files}
)

INCLUDE_DIRECTORIES(
	.
	..
)


IF (PARALLEL_USE_MPI)
	INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF (PARALLEL_USE_MPI)

target_link_libraries (FemLib
	BaseLib
	MathLib
	GeoLib
	MeshLib
	DiscreteLib
)

SET_TARGET_PROPERTIES(FemLib PROPERTIES LINKER_LANGUAGE CXX)
