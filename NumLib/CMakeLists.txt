#Source files grouped by a directory
FILE(GLOB NumLib_Files *.h *.cpp)
SOURCE_GROUP( NumLib FILES ${NumLib_Files} )

ADD_LIBRARY( NumLib STATIC ${NumLib_Files} )

INCLUDE_DIRECTORIES(
	.
	..
)

target_link_libraries (
	NumLib
	BASE
        MathLib
        MeshLib
        FemLib
)
