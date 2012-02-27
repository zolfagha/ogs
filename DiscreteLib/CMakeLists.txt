#Source files grouped by a directory
FILE(GLOB DiscreteLib_Files *.h *.cpp)
SOURCE_GROUP( DiscreteLib FILES ${DiscreteLib_Files} )
FILE(GLOB DiscreteLib_OGS5_Files ogs5/*.h ogs5/*.cpp)
SOURCE_GROUP( DiscreteLib\\OGS5 FILES ${DiscreteLib_OGS5_Files} )


ADD_LIBRARY( DiscreteLib STATIC 
	${DiscreteLib_Files} 
	${DiscreteLib_OGS5_Files} 
)

INCLUDE_DIRECTORIES(
	.
	..
)

target_link_libraries (
	DiscreteLib
	Base
        MathLib
        MeshLib
)
