# Source files
FILE(GLOB GeoLib_Files GeoLib/*.h GeoLib/*.cpp)
SOURCE_GROUP( ogsGeometry FILES ${GeoLib_Files})

FILE(GLOB GeoLib_Core_Files GeoLib/Core/*.h GeoLib/Core/*.cpp)
SOURCE_GROUP( ogsGeometry\\Core FILES ${GeoLib_Core_Files})

FILE(GLOB GeoLib_Shape_Files GeoLib/Shape/*.h GeoLib/Shape/*.cpp)
SOURCE_GROUP( ogsGeometry\\Shape FILES ${GeoLib_Shape_Files})


LIST(APPEND GeoLib_Files
	${GeoLib_Core_Files}
	${GeoLib_Shape_Files}
)
#SET( GeoLib_Files ${GeoLib_Files} PARENT_SCOPE )
#MESSAGE (STATUS "GeoLib files are ${GeoLib_Files}")


# Create the library
#ADD_LIBRARY(ogsGeometry STATIC 
#	${GeoLib_Files} 
#	${GeoLib_Core_Files}
#	${GeoLib_Shape_Files}
#)

#include_directories(
#	.
#	../ogsCore
#)


#target_link_libraries (
#	GeoLib
#	BaseLib
#	MathLib
#)

#SET_TARGET_PROPERTIES(ogsGeometry PROPERTIES LINKER_LANGUAGE CXX)

