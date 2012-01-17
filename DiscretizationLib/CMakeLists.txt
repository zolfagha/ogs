#Source files grouped by a directory
SET ( DiscretizationLib_Files
	DoF.h
)
SOURCE_GROUP( DiscretizationLib FILES ${DiscretizationLib_Files} )

ADD_LIBRARY( DiscretizationLib STATIC ${DiscretizationLib_Files} )

INCLUDE_DIRECTORIES(
	.
	../Base
	../MathLib
)

target_link_libraries (
	DiscretizationLib
	BASE
        MathLib
)
