# Add subdirectories with the projects
INCLUDE( MeshLib/CMakeLists.txt )

#MESSAGE (STATUS "MeshLib files are ${MeshLib_Files}")

ADD_LIBRARY( ogsMesh STATIC 
	${MeshLib_Files}
)

SET_TARGET_PROPERTIES(ogsMesh PROPERTIES LINKER_LANGUAGE CXX)

INCLUDE_DIRECTORIES(
        .
        ../ogsCore
        ../ogsGeometry
        ../extra/Eigen
)

