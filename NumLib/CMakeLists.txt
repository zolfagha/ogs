#Source files grouped by a directory
FILE(GLOB NumLib_Files *.h *.cpp)
SOURCE_GROUP( NumLib FILES ${NumLib_Files} )
FILE(GLOB NumLib_Core_Files Core/*.h Core/*.cpp)
SOURCE_GROUP( NumLib\\Core FILES ${NumLib_Core_Files} )
FILE(GLOB NumLib_TimeStepping_Files TimeStepping/*.h TimeStepping/*.cpp)
SOURCE_GROUP( NumLib\\TimeStepping FILES ${NumLib_TimeStepping_Files} )
FILE(GLOB NumLib_Coupling_Files Coupling/*.h Coupling/*.cpp)
SOURCE_GROUP( NumLib\\Coupling FILES ${NumLib_Coupling_Files} )
FILE(GLOB NumLib_Output_Files Output/*.h Output/*.cpp)
SOURCE_GROUP( NumLib\\Output FILES ${NumLib_Output_Files} )
FILE(GLOB NumLib_Discrete_Files Discrete/*.h Discrete/*.cpp)
SOURCE_GROUP( NumLib\\Discrete FILES ${NumLib_Discrete_Files} )


ADD_LIBRARY( NumLib STATIC 
	${NumLib_Files} 
	${NumLib_Core_Files} 
	${NumLib_TimeStepping_Files} 
	${NumLib_Coupling_Files} 
	${NumLib_Output_Files} 
	${NumLib_Discrete_Files} 
)

INCLUDE_DIRECTORIES(
	.
	..
)

IF (USE_MPI)
	INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF (USE_MPI)

target_link_libraries (
	NumLib
	Base
        MathLib
        MeshLib
	DiscreteLib
)
