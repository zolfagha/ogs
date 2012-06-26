# Source files
FILE(GLOB CUI_Files *.h *.cpp)

#SET( BaseLib_Files ${BaseLib_Files} PARENT_SCOPE )
SOURCE_GROUP( ogsCUI FILES ${CUI_Files} )

#MESSAGE (STATUS "BaseLib files are  ${BaseLib_Files}")

INCLUDE_DIRECTORIES(
        .
        ../extra
)


## Create the library
ADD_EXECUTABLE( ogs6 ${CUI_Files})
TARGET_LINK_LIBRARIES(ogs6 logog)

SET_TARGET_PROPERTIES(ogs6 PROPERTIES LINKER_LANGUAGE CXX)

