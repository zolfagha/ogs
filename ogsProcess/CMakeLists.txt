# Add subdirectories with the projects
INCLUDE( ProcessLib/CMakeLists.txt )
INCLUDE( ApplicationLib/CMakeLists.txt )

ADD_LIBRARY( ogsProcess STATIC 
	${ProcessLib_Files}
    ${ApplicationLib_Files}
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

