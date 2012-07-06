## pthread ##
SET ( CMAKE_THREAD_PREFER_PTHREAD On )
FIND_PACKAGE( Threads )
IF ( CMAKE_USE_PTHREADS_INIT )
        SET (HAVE_PTHREADS TRUE)
        MESSAGE (STATUS "pthread library found." )
ENDIF (CMAKE_USE_PTHREADS_INIT )


# Source files
FILE(GLOB CUI_Files *.h *.cpp)
SOURCE_GROUP( ogsCUI FILES ${CUI_Files} )


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
        ../ogsFileIO
        ../ogsSimulator
)


## Create the library
ADD_EXECUTABLE( ogs6 ${CUI_Files})
TARGET_LINK_LIBRARIES( ogs6 
    ogsFileIO
    ogsProcess
    ogsSimulator
    logog
    ${CMAKE_THREAD_LIBS_INIT}
    )

IF(USE_LIS)
    INCLUDE_DIRECTORIES( ${LIS_INCLUDE_DIR} )
ENDIF(USE_LIS)

IF (USE_LIS)
    TARGET_LINK_LIBRARIES( ogs6 ${LIS_LIBRARIES} )
ENDIF (USE_LIS)

IF (USE_BLAS_LAPACK)
    TARGET_LINK_LIBRARIES( ogs6 ${BLAS_LIBRARIES} )
    TARGET_LINK_LIBRARIES( ogs6 ${LAPACK_LIBRARIES} )
ENDIF (USE_BLAS_LAPACK)

IF(USE_PETSC)
        TARGET_LINK_LIBRARIES( ogs6  ${PETSC_LIBRARIES} )
    IF (MSVC)
        SET_TARGET_PROPERTIES(ogs6 PROPERTIES
                LINK_FLAGS "
            /NODEFAULTLIB:\"libcmt.lib\"
                ")
    ENDIF(MSVC)
ENDIF(USE_PETSC)

IF (USE_MPI)
    INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
        TARGET_LINK_LIBRARIES( ogs6 ${MPI_C_LIBRARIES} )
ENDIF (USE_MPI)


SET_TARGET_PROPERTIES(ogs6 PROPERTIES LINKER_LANGUAGE CXX)

