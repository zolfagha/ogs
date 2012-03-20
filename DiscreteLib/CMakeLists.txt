#Source files grouped by a directory
# ./
FILE(GLOB DiscreteLib_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h *.cpp)
SOURCE_GROUP( DiscreteLib FILES ${DiscreteLib_Files} )
# ./Core
FILE(GLOB DiscreteLib_Core_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} Core/*.h Core/*.cpp)
SOURCE_GROUP( DiscreteLib FILES ${DiscreteLib_Core_Files} )
# ./Linear equation
FILE(GLOB DiscreteLib_LinearEquation_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} LinearEquation/*.h LinearEquation/*.cpp)
SOURCE_GROUP( DiscreteLib FILES ${DiscreteLib_LinearEquation_Files} )
# ./EquationId
FILE(GLOB DiscreteLib_EquationId_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} EquationId/*.h EquationId/*.cpp)
SOURCE_GROUP( DiscreteLib FILES ${DiscreteLib_EquationId_Files} )
# ./OGS5
IF(USE_MPI)
	FILE(GLOB DiscreteLib_OGS5_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ogs5/*.h ogs5/*.cpp)
	SOURCE_GROUP( DiscreteLib\\OGS5 FILES ${DiscreteLib_OGS5_Files} )
ENDIF()
# ./Utils
FILE(GLOB DiscreteLib_Utils_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} Utils/*.h Utils/*.cpp)
SOURCE_GROUP( DiscreteLib FILES ${DiscreteLib_Utils_Files} )
# ./Assembler
FILE(GLOB DiscreteLib_Assembler_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} Assembler/*.h Assembler/*.cpp)
SOURCE_GROUP( DiscreteLib FILES ${DiscreteLib_Assembler_Files} )
# ./DDC
FILE(GLOB DiscreteLib_DDC_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} DDC/*.h DDC/*.cpp)
SOURCE_GROUP( DiscreteLib FILES ${DiscreteLib_DDC_Files} )
IF(NOT USE_LIS OR NOT USE_MPI)
    LIST(REMOVE_ITEM DiscreteLib_DDC_Files 
        DDC/LisMPIDiscreteVector.h 
        DDC/LisMPIDiscreteVector.cpp
        DDC/LisMPILinearSystem.h 
        DDC/LisDiscreteSystem.h
        DDC/LisDiscreteSystem.cpp
        )
ENDIF()
# ./OpenMP
FILE(GLOB DiscreteLib_OpenMP_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} OpenMP/*.h OpenMP/*.cpp)
SOURCE_GROUP( DiscreteLib FILES ${DiscreteLib_OpenMP_Files} )

ADD_LIBRARY( DiscreteLib STATIC 
	${DiscreteLib_Files} 
    ${DiscreteLib_Core_Files} 
    ${DiscreteLib_LinearEquation_Files}
    ${DiscreteLib_EquationId_Files} 
    ${DiscreteLib_Assembler_Files} 
    ${DiscreteLib_Utils_Files} 
    ${DiscreteLib_DDC_Files} 
    ${DiscreteLib_OpenMP_Files} 
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
