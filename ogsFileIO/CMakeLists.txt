#
INCLUDE(ogs5/CMakeLists.txt)
INCLUDE(MeshIO/CMakeLists.txt)
INCLUDE(OutputIO/CMakeLists.txt)

# Create the library
ADD_LIBRARY(ogsFileIO STATIC 
	${OGS5_FILES}
	${MESHIO_FILES}
	${OUTPUTIO_FILES}
)

include_directories(
	.
)
