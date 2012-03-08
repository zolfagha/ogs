#Source files grouped by a directory
#SET ( FdmLib_Files
#	FDM.h
#)
FILE(GLOB FdmLib_Files *.h *.cpp)
SOURCE_GROUP( FdmLib FILES ${FdmLib_Files} )

ADD_LIBRARY( FdmLib STATIC ${FdmLib_Files} )

INCLUDE_DIRECTORIES(
	.
	../Base
	../MathLib
)

target_link_libraries (
	FdmLib
	Base
        MathLib
)
SET_TARGET_PROPERTIES(Base PROPERTIES LINKER_LANGUAGE CXX)
