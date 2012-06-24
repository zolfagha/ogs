# Source files
FILE(GLOB Base_Files *.h *.cpp)
SOURCE_GROUP( BaseLib FILES ${Base_Files} )

# Create the library
ADD_LIBRARY( BaseLib STATIC ${Base_Files})

SET_TARGET_PROPERTIES(BaseLib PROPERTIES LINKER_LANGUAGE CXX)

INCLUDE_DIRECTORIES(
        .
)

# Add logog subdirectory and group its targets in a Visual Studio folder
ADD_SUBDIRECTORY(logog)
SET_PROPERTY(TARGET Continuous PROPERTY FOLDER "logog")
SET_PROPERTY(TARGET logog PROPERTY FOLDER "logog")
IF (DOXYGEN_FOUND)
	SET_PROPERTY(TARGET logog-doc PROPERTY FOLDER "logog")
ENDIF ()
SET_PROPERTY(TARGET test-logog PROPERTY FOLDER "logog")
