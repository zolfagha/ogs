
FILE(GLOB Test_Files *.h *.cpp)
SOURCE_GROUP( Test FILES ${Test_Files} )
FILE(GLOB Test_Solution1_Files Solution1/*.h Solution1/*.cpp)
SOURCE_GROUP( Test\\Solution1 FILES ${Test_Solution1_Files} )
FILE(GLOB Test_Geo_Material_Files Geo/Material/*.h Geo/Material/*.cpp)
SOURCE_GROUP( Test\\Geo\\Material FILES ${Test_Geo_Material_Files} )
FILE(GLOB Test_Geo_Equation_Files Geo/Equation/*.h Geo/Equation/*.cpp)
SOURCE_GROUP( Test\\Geo\\Equation FILES ${Test_Geo_Equation_Files} )
FILE(GLOB Test_Geo_Model_Files Geo/Model/*.h Geo/Model/*.cpp)
SOURCE_GROUP( Test\\Geo\\Model FILES ${Test_Geo_Model_Files} )
FILE(GLOB Test_ExactSolution_Files ExactSolution/*.h ExactSolution/*.cpp)
SOURCE_GROUP( Test\\ExactSolution FILES ${Test_ExactSolution_Files} )
FILE(GLOB Test_FDM_Files FDM/*.h FDM/*.cpp)
SOURCE_GROUP( Test\\FDM FILES ${Test_FDM_Files} )

INCLUDE_DIRECTORIES(
	.
	..
)

# Create the executable
ADD_EXECUTABLE( tests
        ${Test_Files}
        ${Test_Solution1_Files}
        ${Test_Geo_Material_Files}
        ${Test_Geo_Equation_Files}
        ${Test_Geo_Model_Files}
        ${Test_ExactSolution_Files}
        ${Test_FDM_Files}
)

TARGET_LINK_LIBRARIES ( tests
	Base
	MathLib
	GeoLib
	MeshLib
	DiscreteLib
	#FdmLib
	FemLib
	NumLib
)


# Link external libraries
TARGET_LINK_LIBRARIES( tests ${GTEST_BOTH_LIBRARIES} )

IF (WIN32)
        TARGET_LINK_LIBRARIES( tests Winmm.lib )
ENDIF (WIN32)

IF (LIS_FOUND)
	TARGET_LINK_LIBRARIES( tests ${LIS_LIBRARIES} )
ENDIF (LIS_FOUND)

IF (USE_BLAS_LAPACK)
	TARGET_LINK_LIBRARIES( tests ${BLAS_LIBRARIES} )
	TARGET_LINK_LIBRARIES( tests ${LAPACK_LIBRARIES} )
ENDIF (USE_BLAS_LAPACK)

IF(USE_PETSC)
        TARGET_LINK_LIBRARIES( tests  ${PETSC_LIBRARIES} )
	IF (MSVC)
		SET_TARGET_PROPERTIES(tests PROPERTIES
		    	LINK_FLAGS "
			/NODEFAULTLIB:\"libcmt.lib\"
		    	")
	ENDIF( MSVC)
ENDIF(USE_PETSC)

IF (USE_MPI)
	INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
        TARGET_LINK_LIBRARIES( tests ${MPI_C_LIBRARIES} )
ENDIF (USE_MPI)


SET_TARGET_PROPERTIES(tests PROPERTIES LINKER_LANGUAGE CXX)

