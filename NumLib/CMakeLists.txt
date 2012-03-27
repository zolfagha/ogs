#Source files grouped by a directory
FILE(GLOB NumLib_Files *.h *.cpp)
SOURCE_GROUP( NumLib FILES ${NumLib_Files} )
FILE(GLOB NumLib_Core_Files Core/*.h Core/*.cpp)
SOURCE_GROUP( NumLib\\Core FILES ${NumLib_Core_Files} )
FILE(GLOB NumLib_TimeStepping_Files TimeStepping/*.h TimeStepping/*.cpp)
SOURCE_GROUP( NumLib\\TimeStepping FILES ${NumLib_TimeStepping_Files} )
FILE(GLOB NumLib_Coupling_Files Coupling/*.h Coupling/*.cpp)
SOURCE_GROUP( NumLib\\Coupling FILES ${NumLib_Coupling_Files} )
FILE(GLOB NumLib_Coupling_Algorithm_Files Coupling/Algorithm/*.h TransientCoupling/Algorithm/*.cpp)
SOURCE_GROUP( NumLib\\Coupling\\Algorithm FILES ${NumLib_Coupling_Algorithm_Files} )
FILE(GLOB NumLib_TransientCoupling_Files TransientCoupling/*.h TransientCoupling/*.cpp)
SOURCE_GROUP( NumLib\\TransientCoupling FILES ${NumLib_TransientCoupling_Files} )
FILE(GLOB NumLib_Output_Files Output/*.h Output/*.cpp)
SOURCE_GROUP( NumLib\\Output FILES ${NumLib_Output_Files} )
FILE(GLOB NumLib_TransientAssembler_Files TransientAssembler/*.h TransientAssembler/*.cpp)
SOURCE_GROUP( NumLib\\TransientAssembler FILES ${NumLib_TransientAssembler_Files} )



ADD_LIBRARY( NumLib STATIC 
	${NumLib_Files} 
	${NumLib_Core_Files} 
	${NumLib_TimeStepping_Files} 
	${NumLib_Coupling_Files} 
    ${NumLib_Coupling_Algorithm_Files} 
    ${NumLib_TransientCoupling_Files} 
	${NumLib_Output_Files} 
	${NumLib_TransientAssembler_Files} 
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
SET_TARGET_PROPERTIES(NumLib PROPERTIES LINKER_LANGUAGE CXX)
