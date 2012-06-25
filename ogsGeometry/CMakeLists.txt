# Add subdirectories with the projects
#ADD_SUBDIRECTORY( GeoLib)
#SOURCE_GROUP( GeoLib FILES ${GeoLib_Files} )
INCLUDE( GeoLib/CMakeLists.txt )

#MESSAGE (STATUS "GeoLib files are ${GeoLib_Files}")

ADD_LIBRARY( ogsGeometry STATIC 
	${GeoLib_Files}
)

SET_TARGET_PROPERTIES(ogsGeometry PROPERTIES LINKER_LANGUAGE CXX)

INCLUDE_DIRECTORIES(
        .
        ../ogsCore
        ../extra/Eigen
)

