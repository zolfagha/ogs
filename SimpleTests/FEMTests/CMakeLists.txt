
INCLUDE_DIRECTORIES(
        .
	../../Base/
	../../MathLib/
	../../GeoLib/
	../../MeshLib/
	../../FemLib/
)

# Create the executable
ADD_EXECUTABLE( FEMTests
        FEMTests.cpp
        ${SOURCES}
        ${HEADERS}
)

#IF (WIN32)
#        TARGET_LINK_LIBRARIES(MatMult Winmm.lib)
#ENDIF (WIN32)

TARGET_LINK_LIBRARIES ( FEMTests
	Base
	MathLib
	GeoLib
#	MeshLib
#	FemLib
)

