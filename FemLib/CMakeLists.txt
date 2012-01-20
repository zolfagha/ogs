#Source files grouped by a directory
#SET ( FemLib_Files
#	Fem.h
#	FemFunction.h
#	MeshSparseTable.h
#	Integration.h
#	Mapping.h
#	WeightFunction.h
#	Test.cpp
#
FILE(GLOB FemLib_Files *.h *.cpp)
SOURCE_GROUP( FemLib FILES ${FemLib_Files})
FILE(GLOB FemLib_Element_Files Element/*.h *.cpp)
SOURCE_GROUP( FemLib\\Element FILES ${FemLib_Element_Files})
FILE(GLOB FemLib_BC_Files BC/*.h *.cpp)
SOURCE_GROUP( FemLib\\BC FILES ${FemLib_BC_Files})

ADD_LIBRARY(FemLib STATIC
	${FemLib_Files}
	${FemLib_Element_Files}
	${FemLib_BC_Files}
)

INCLUDE_DIRECTORIES(
	.
	..
#	../Base
#	../GeoLib
#	../MathLib
#	../MeshLib
)


IF (PARALLEL_USE_MPI)
	INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF (PARALLEL_USE_MPI)

target_link_libraries (
	FemLib
	Base
        MathLib
	GeoLib
	MeshLib
)
