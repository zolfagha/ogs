# Add subdirectories with the projects
INCLUDE( ProcessLib/CMakeLists.txt )

ADD_LIBRARY( ogsProcess STATIC 
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

