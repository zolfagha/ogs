# Source files
FILE(GLOB ProcessLib_Files ProcessLib/*.h ProcessLib/*.cpp)

#SET( BaseLib_Files ${BaseLib_Files} PARENT_SCOPE )
SOURCE_GROUP( ogsProcess\\ProcessLib FILES ${ProcessLib_Files} )

#MESSAGE (STATUS "BaseLib files are  ${BaseLib_Files}")

## Create the library
#ADD_LIBRARY( BaseLib STATIC ${Base_Files})

#SET_TARGET_PROPERTIES(BaseLib PROPERTIES LINKER_LANGUAGE CXX)

#INCLUDE_DIRECTORIES(
#        .
#)

