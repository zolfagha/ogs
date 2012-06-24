#Source files grouped by a directory
FILE(GLOB NumLib_Files *.h *.cpp)
SOURCE_GROUP( NumLib FILES ${NumLib_Files} )
FILE(GLOB NumLib_TimeStepping_Files TimeStepping/*.h TimeStepping/*.cpp)
SOURCE_GROUP( NumLib\\TimeStepping FILES ${NumLib_TimeStepping_Files} )
FILE(GLOB NumLib_TransientCoupling_Files TransientCoupling/*.h TransientCoupling/*.cpp)
SOURCE_GROUP( NumLib\\TransientCoupling FILES ${NumLib_TransientCoupling_Files} )
FILE(GLOB NumLib_FileIO_Files FileIO/*.h FileIO/*.cpp)
SOURCE_GROUP( NumLib\\FileIO FILES ${NumLib_FileIO_Files} )
FILE(GLOB NumLib_TransientAssembler_Files TransientAssembler/*.h TransientAssembler/*.cpp)
SOURCE_GROUP( NumLib\\TransientAssembler FILES ${NumLib_TransientAssembler_Files} )
FILE(GLOB NumLib_Nonlinear_Files Nonlinear/*.h Nonlinear/*.cpp)
SOURCE_GROUP( NumLib\\Nonlinear FILES ${NumLib_Nonlinear_Files} )
#Function/
FILE(GLOB NumLib_Function_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
	Function/*.h 
	Function/*.cpp)
SOURCE_GROUP( NumLib\\Function FILES ${NumLib_Function_Files})
#IOSystem/
FILE(GLOB NumLib_IOSystem_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
    IOSystem/*.h 
    IOSystem/*.cpp)
SOURCE_GROUP( NumLib\\IOSystem FILES ${NumLib_IOSystem_Files})
#Coupling/
FILE(GLOB NumLib_Coupling_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
    Coupling/*.h 
    Coupling/*.cpp)
SOURCE_GROUP( NumLib\\Coupling FILES ${NumLib_Coupling_Files})
#Coupling/Algorithm
FILE(GLOB NumLib_Coupling_Algorithm_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
    Coupling/Algorithm/*.h 
    Coupling/Algorithm/*.cpp)
SOURCE_GROUP( NumLib\\Coupling\\Algorithm FILES ${NumLib_Coupling_Algorithm_Files})
#SystemOfEquations/
FILE(GLOB NumLib_SystemOfEquations_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
    SystemOfEquations/*.h 
    SystemOfEquations/*.cpp)
SOURCE_GROUP( NumLib\\SystemOfEquations FILES ${NumLib_SystemOfEquations_Files})

ADD_LIBRARY( NumLib STATIC 
	${NumLib_Files} 
#	${NumLib_Core_Files} 
    ${NumLib_Nonlinear_Files} 
    ${NumLib_FileIO_Files} 
	${NumLib_TimeStepping_Files} 
    ${NumLib_TransientCoupling_Files} 
	${NumLib_TransientAssembler_Files} 
    ${NumLib_Coupling_Files}
    ${NumLib_Coupling_Algorithm_Files}
    ${NumLib_Function_Files}
    ${NumLib_IOSystem_Files}
    ${NumLib_SystemOfEquations_Files}
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
	BaseLib
        NumLib
        MeshLib
	DiscreteLib
)
SET_TARGET_PROPERTIES(NumLib PROPERTIES LINKER_LANGUAGE CXX)
