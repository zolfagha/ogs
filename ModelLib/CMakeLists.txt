#Source files grouped by a directory
FILE(GLOB ModelLib_Files *.h *.cpp)
SOURCE_GROUP( ModelLib FILES ${ModelLib_Files})

ADD_LIBRARY(ModelLib STATIC ${ModelLib_Files})

INCLUDE_DIRECTORIES(
	.
	..
)


target_link_libraries (
	ModelLib
	NumLib
)
