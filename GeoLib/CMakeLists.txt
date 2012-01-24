# Source files
FILE(GLOB GeoLib_Files *.h *.cpp)
SOURCE_GROUP( GeoLib FILES ${GeoLib_Files})

FILE(GLOB GeoLib_Core_Files Core/*.h Core/*.cpp)
SOURCE_GROUP( GeoLib\\Core FILES ${GeoLib_Core_Files})

FILE(GLOB GeoLib_Shape_Files Shape/*.h Shape/*.cpp)
SOURCE_GROUP( GeoLib\\Shape FILES ${GeoLib_Shape_Files})

# Create the library
ADD_LIBRARY(GeoLib STATIC 
	${GeoLib_Files} 
	${GeoLib_Core_Files}
	${GeoLib_Shape_Files}
)

include_directories(
	.
	..
)


target_link_libraries (
	GeoLib
	Base
	MathLib
)

