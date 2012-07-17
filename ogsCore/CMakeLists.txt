# Add subdirectories with the projects
#ADD_SUBDIRECTORY( BaseLib )
#SOURCE_GROUP( BaseLib FILES ${BaseLib_Files} )
INCLUDE( BaseLib/CMakeLists.txt )
#ADD_SUBDIRECTORY( GeoLib)
#SOURCE_GROUP( GeoLib FILES ${GeoLib_Files} )
#INCLUDE( GeoLib/CMakeLists.txt )
#ADD_SUBDIRECTORY( MathLib )
#SOURCE_GROUP( MathLib FILES ${MathLib_Files} )
INCLUDE( MathLib/CMakeLists.txt )

#MESSAGE (STATUS "GeoLib files are ${GeoLib_Files}")
#MESSAGE (STATUS "MathLib files are ${MathLib_Files}")

ADD_LIBRARY( ogsCore STATIC 
	${BaseLib_Files}
#	${GeoLib_Files}
	${MathLib_Files}
)

SET_TARGET_PROPERTIES(ogsCore PROPERTIES LINKER_LANGUAGE CXX)

INCLUDE_DIRECTORIES(
        .
        ../extra/Eigen
        ../extra/tinyxml2
)

