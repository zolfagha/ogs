
# Source files
FILE(GLOB Simulator_Files *.h *.cpp)

#SET( BaseLib_Files ${BaseLib_Files} PARENT_SCOPE )
SOURCE_GROUP( ogsOGS6 FILES ${Simulator_Files} )

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

ADD_LIBRARY( ogsOGS6 STATIC 
	${Simulator_Files}
)

target_link_libraries( ogsOGS6
    ogsFileIO
)

SET_TARGET_PROPERTIES(ogsOGS6 PROPERTIES LINKER_LANGUAGE CXX)

