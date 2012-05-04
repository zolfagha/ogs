#Source files grouped by a directory
FILE(GLOB NumLib_Files *.h *.cpp)
SOURCE_GROUP( NumLib FILES ${NumLib_Files} )
FILE(GLOB NumLib_TimeStepping_Files TimeStepping/*.h TimeStepping/*.cpp)
SOURCE_GROUP( NumLib\\TimeStepping FILES ${NumLib_TimeStepping_Files} )
FILE(GLOB NumLib_TransientCoupling_Files TransientCoupling/*.h TransientCoupling/*.cpp)
SOURCE_GROUP( NumLib\\TransientCoupling FILES ${NumLib_TransientCoupling_Files} )
FILE(GLOB NumLib_Output_Files Output/*.h Output/*.cpp)
SOURCE_GROUP( NumLib\\Output FILES ${NumLib_Output_Files} )
FILE(GLOB NumLib_TransientAssembler_Files TransientAssembler/*.h TransientAssembler/*.cpp)
SOURCE_GROUP( NumLib\\TransientAssembler FILES ${NumLib_TransientAssembler_Files} )
FILE(GLOB NumLib_Nonlinear_Files Nonlinear/*.h Nonlinear/*.cpp)
SOURCE_GROUP( NumLib\\Nonlinear FILES ${NumLib_Nonlinear_Files} )
FILE(GLOB NumLib_TXFunction_Files TXFunction/*.h TXFunction/*.cpp)
SOURCE_GROUP( NumLib\\TXFunction FILES ${NumLib_TXFunction_Files} )


ADD_LIBRARY( NumLib STATIC 
	${NumLib_Files} 
#	${NumLib_Core_Files} 
    ${NumLib_Nonlinear_Files} 
    ${NumLib_Output_Files} 
	${NumLib_TimeStepping_Files} 
    ${NumLib_TransientCoupling_Files} 
	${NumLib_TransientAssembler_Files} 
    ${NumLib_TXFunction_Files} 
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
