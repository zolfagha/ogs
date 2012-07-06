
# Source files
FILE(GLOB Simulator_Files *.h *.cpp)

#SET( BaseLib_Files ${BaseLib_Files} PARENT_SCOPE )
SOURCE_GROUP( ogsSimulator FILES ${Simulator_Files} )

INCLUDE_DIRECTORIES(
        .
        ../ogsCore
        ../ogsGeometry
        ../ogsMesh
        ../ogsNumerics
        ../ogsPhysics
        ../ogsFileIO
        ../ogsProcess
        ../extra/Eigen
)

ADD_LIBRARY( ogsSimulator STATIC 
	${Simulator_Files}
)

target_link_libraries( ogsSimulator
    ogsGeometry
    ogsMesh
    ogsNumerics
    ogsFileIO
    ogsProcess
)

SET_TARGET_PROPERTIES(ogsSimulator PROPERTIES LINKER_LANGUAGE CXX)

