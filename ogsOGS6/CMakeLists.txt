
# Source files
FILE(GLOB Simulator_Files *.h *.cpp)

#SET( BaseLib_Files ${BaseLib_Files} PARENT_SCOPE )
SOURCE_GROUP( ogsOGS6 FILES ${Simulator_Files} )

ADD_LIBRARY( ogsOGS6 STATIC 
	${Simulator_Files}
)

SET_TARGET_PROPERTIES(ogsOGS6 PROPERTIES LINKER_LANGUAGE CXX)

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

