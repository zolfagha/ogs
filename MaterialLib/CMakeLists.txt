#Source files grouped by a directory
SET ( MaterialLib_Files
	PorousMedium.h
	SolidMechanics.h
)
SOURCE_GROUP( MaterialLib FILES ${MaterialLib_Files} )

ADD_LIBRARY( MaterialLib STATIC ${MaterialLib_Files} )

INCLUDE_DIRECTORIES(
	.
	..
)

target_link_libraries (
	MaterialLib
	Base
        MathLib
)
