# Add subdirectories with the projects
INCLUDE( ProcessLib/CMakeLists.txt )

FILE(GLOB Util_Files *.h *.cpp)
SOURCE_GROUP( ogsProcess FILES ${Util_Files} )

ADD_LIBRARY( ogsProcess STATIC 
    ${Util_Files}
	${ProcessLib_Files}
)

SET_TARGET_PROPERTIES(ogsProcess PROPERTIES LINKER_LANGUAGE CXX)

INCLUDE_DIRECTORIES(
        .
        ../ogsCore
        ../ogsGeometry
        ../ogsMesh
        ../ogsNumerics
        ../ogsPhysics
        ../extra/Eigen
)
