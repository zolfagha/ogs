## pthread ##
SET ( CMAKE_THREAD_PREFER_PTHREAD On )
FIND_PACKAGE( Threads )
IF ( CMAKE_USE_PTHREADS_INIT )
        SET (HAVE_PTHREADS TRUE)
        MESSAGE (STATUS "pthread library found." )
ENDIF (CMAKE_USE_PTHREADS_INIT )


# Source files
FILE(GLOB CUI_Files *.h *.cpp)

#SET( BaseLib_Files ${BaseLib_Files} PARENT_SCOPE )
SOURCE_GROUP( ogsCUI FILES ${CUI_Files} )

#MESSAGE (STATUS "BaseLib files are  ${BaseLib_Files}")

INCLUDE_DIRECTORIES(
        .
        ../extra
        ../extra/Eigen
        ../ogsCore
        ../ogsGeometry
        ../ogsMesh
        ../ogsNumerics
        ../ogsPhysics
        ../ogsProcess
)


## Create the library
ADD_EXECUTABLE( ogs6 ${CUI_Files})
TARGET_LINK_LIBRARIES(ogs6 
    ogsProcess
    logog
    ${CMAKE_THREAD_LIBS_INIT}
    )

SET_TARGET_PROPERTIES(ogs6 PROPERTIES LINKER_LANGUAGE CXX)

