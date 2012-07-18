#
INCLUDE(CommonIO/CMakeLists.txt)
INCLUDE(GeoIO/CMakeLists.txt)
INCLUDE(MeshIO/CMakeLists.txt)
INCLUDE(FemIO/CMakeLists.txt)
INCLUDE(OutputIO/CMakeLists.txt)

# Create the library
ADD_LIBRARY(ogsFileIO STATIC 
    ${COMMONIO_FILES}
	${GEOIO_FILES}
	${MESHIO_FILES}
    ${FEMIO_FILES}
	${OUTPUTIO_FILES}
)

include_directories(
	.
    ../extra/Eigen
	../ogsCore
    ../ogsGeometry
	../ogsMesh
    ../ogsNumerics
)

target_link_libraries( ogsFileIO
    ogsCore
    ogsGeometry
    ogsMesh
    ogsNumerics
)

SET_TARGET_PROPERTIES(ogsFileIO PROPERTIES LINKER_LANGUAGE CXX)
