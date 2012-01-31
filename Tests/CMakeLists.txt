
FILE(GLOB Test_Files *.h *.cpp)

INCLUDE_DIRECTORIES(
	.
	..
)

# Create the executable
ADD_EXECUTABLE( tests
        ${Test_Files}
)

TARGET_LINK_LIBRARIES(tests ${GTEST_BOTH_LIBRARIES} )


TARGET_LINK_LIBRARIES ( tests
	Base
	MathLib
	GeoLib
	MeshLib
#	FdmLib
	FemLib
	NumLib
)

