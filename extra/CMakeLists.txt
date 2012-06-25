# Add logog subdirectory and group its targets in a Visual Studio folder
ADD_SUBDIRECTORY(logog)
SET_PROPERTY(TARGET Continuous PROPERTY FOLDER "logog")
SET_PROPERTY(TARGET logog PROPERTY FOLDER "logog")
IF (DOXYGEN_FOUND)
	SET_PROPERTY(TARGET logog-doc PROPERTY FOLDER "logog")
ENDIF ()
SET_PROPERTY(TARGET test-logog PROPERTY FOLDER "logog")
