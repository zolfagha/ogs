# Source files
FILE(GLOB Base_Files *.h *.cpp)
SOURCE_GROUP( Base FILES ${Base_Files} )

# Create the library
ADD_LIBRARY( Base STATIC ${Base_Files})

SET_TARGET_PROPERTIES(Base PROPERTIES LINKER_LANGUAGE CXX)

INCLUDE_DIRECTORIES(
        .
)

