# Add subdirectories with the projects
INCLUDE( MaterialLib/CMakeLists.txt )
INCLUDE( ProcessLib/CMakeLists.txt )

ADD_LIBRARY( ogsPhysics STATIC 
	${MaterialLib_Files}
	${ProcessLib_Files}
)

SET_TARGET_PROPERTIES(ogsPhysics PROPERTIES LINKER_LANGUAGE CXX)

INCLUDE_DIRECTORIES(
        .
        ../ogsCore
        ../ogsGeometry
        ../ogsMesh
        ../ogsNumerics
        ../extra/Eigen
)

