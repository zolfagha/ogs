#Source files grouped by a directory
# ./
FILE(GLOB ProcessLib_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ProcessLib/*.h ProcessLib/*.cpp)
SOURCE_GROUP( ogsPhysics\\ProcessLib FILES ${ProcessLib_Files} )
# ./Core
#FILE(GLOB DiscreteLib_Core_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} DiscreteLib/Core/*.h DiscreteLib/Core/*.cpp)
#SOURCE_GROUP( ogsNumerics\\DiscreteLib\\Core FILES ${DiscreteLib_Core_Files} )


#
#LIST(APPEND MaterialLib_Files
#    ${MaterialLib_Files} 
#) 

