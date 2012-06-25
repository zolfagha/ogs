#Source files grouped by a directory
# ./
FILE(GLOB DiscreteLib_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} DiscreteLib/*.h DiscreteLib/*.cpp)
SOURCE_GROUP( ogsNumerics\\DiscreteLib FILES ${DiscreteLib_Files} )
# ./Core
FILE(GLOB DiscreteLib_Core_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} DiscreteLib/Core/*.h DiscreteLib/Core/*.cpp)
SOURCE_GROUP( ogsNumerics\\DiscreteLib\\Core FILES ${DiscreteLib_Core_Files} )
# ./Linear equation
FILE(GLOB DiscreteLib_LinearEquation_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} DiscreteLib/LinearEquation/*.h DiscreteLib/LinearEquation/*.cpp)
SOURCE_GROUP( ogsNumerics\\DiscreteLib\\LinearEquation FILES ${DiscreteLib_LinearEquation_Files} )
# ./EquationId
FILE(GLOB DiscreteLib_EquationId_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} DiscreteLib/EquationId/*.h DiscreteLib/EquationId/*.cpp)
SOURCE_GROUP( ogsNumerics\\DiscreteLib\\EquationId FILES ${DiscreteLib_EquationId_Files} )
# ./OGS5
IF(USE_MPI)
	FILE(GLOB DiscreteLib_OGS5_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} DiscreteLib/ogs5/*.h DiscreteLib/ogs5/*.cpp)
	SOURCE_GROUP( ogsNumerics\\DiscreteLib\\OGS5 FILES ${DiscreteLib_OGS5_Files} )
ENDIF()
# ./Utils
FILE(GLOB DiscreteLib_Utils_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} DiscreteLib/Utils/*.h DiscreteLib/Utils/*.cpp)
SOURCE_GROUP( ogsNumerics\\DiscreteLib\\Utils FILES ${DiscreteLib_Utils_Files} )
# ./Assembler
FILE(GLOB DiscreteLib_Assembler_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} DiscreteLib/Assembler/*.h DiscreteLib/Assembler/*.cpp)
SOURCE_GROUP( ogsNumerics\\DiscreteLib\\Assembler FILES ${DiscreteLib_Assembler_Files} )
# ./DDC
FILE(GLOB DiscreteLib_DDC_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} DiscreteLib/DDC/*.h DiscreteLib/DDC/*.cpp)
SOURCE_GROUP( ogsNumerics\\DiscreteLib\\DDC FILES ${DiscreteLib_DDC_Files} )
IF(NOT USE_LIS OR NOT USE_MPI)
    LIST(REMOVE_ITEM DiscreteLib_DDC_Files 
        DiscreteLib/DDC/LisMPIDiscreteVector.h 
        DiscreteLib/DDC/LisMPIDiscreteVector.cpp
        DiscreteLib/DDC/LisMPILinearSystem.h 
        DiscreteLib/DDC/LisDiscreteSystem.h
        DiscreteLib/DDC/LisDiscreteSystem.cpp
        )
ENDIF()
# ./OpenMP
FILE(GLOB DiscreteLib_OpenMP_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} DiscreteLib/OpenMP/*.h DiscreteLib/OpenMP/*.cpp)
SOURCE_GROUP( ogsNumerics\\DiscreteLib\\OpenMP FILES ${DiscreteLib_OpenMP_Files} )

#
LIST(APPEND DiscreteLib_Files
    ${DiscreteLib_Core_Files} 
    ${DiscreteLib_LinearEquation_Files}
    ${DiscreteLib_EquationId_Files} 
    ${DiscreteLib_Assembler_Files} 
    ${DiscreteLib_Utils_Files} 
    ${DiscreteLib_DDC_Files} 
    ${DiscreteLib_OpenMP_Files} 
	${DiscreteLib_OGS5_Files} 
)


#ADD_LIBRARY( DiscreteLib STATIC 
#	${DiscreteLib_Files} 
#    ${DiscreteLib_Core_Files} 
#    ${DiscreteLib_LinearEquation_Files}
#    ${DiscreteLib_EquationId_Files} 
#    ${DiscreteLib_Assembler_Files} 
#    ${DiscreteLib_Utils_Files} 
#    ${DiscreteLib_DDC_Files} 
#    ${DiscreteLib_OpenMP_Files} 
#	${DiscreteLib_OGS5_Files} 
#)
#
#INCLUDE_DIRECTORIES(
#	.
#	..
#	../../ogsCore
#)
#
#IF (USE_MPI)
#	INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
#ENDIF (USE_MPI)
#
#target_link_libraries (
#	DiscreteLib
#	BaseLib
#        MathLib
#        MeshLib
#)
#
#SET_TARGET_PROPERTIES(DiscreteLib PROPERTIES LINKER_LANGUAGE CXX)
#