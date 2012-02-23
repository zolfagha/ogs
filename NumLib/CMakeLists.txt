#Source files grouped by a directory
FILE(GLOB NumLib_Files *.h *.cpp)
SOURCE_GROUP( NumLib FILES ${NumLib_Files} )
FILE(GLOB NumLib_Core_Files Core/*.h Core/*.cpp)
SOURCE_GROUP( NumLib\\Core FILES ${NumLib_Core_Files} )
FILE(GLOB NumLib_TimeStepping_Files TimeStepping/*.h TimeStepping/*.cpp)
SOURCE_GROUP( NumLib\\TimeStepping FILES ${NumLib_TimeStepping_Files} )
FILE(GLOB NumLib_Coupling_Files Coupling/*.h Coupling/*.cpp)
SOURCE_GROUP( NumLib\\Coupling FILES ${NumLib_Coupling_Files} )
FILE(GLOB NumLib_Output_Files Output/*.h Output/*.cpp)
SOURCE_GROUP( NumLib\\Output FILES ${NumLib_Output_Files} )
FILE(GLOB NumLib_Discrete_Files Discrete/*.h Discrete/*.cpp)
SOURCE_GROUP( NumLib\\Discrete FILES ${NumLib_Discrete_Files} )
FILE(GLOB NumLib_DomainDecomposition_Files DomainDecomposition/*.h DomainDecomposition/*.cpp)
SOURCE_GROUP( NumLib\\DomainDecomposition FILES ${NumLib_DomainDecomposition_Files} )
FILE(GLOB NumLib_DomainDecomposition_OGS5_Files DomainDecomposition/ogs5/*.h DomainDecomposition/ogs5/*.cpp)
SOURCE_GROUP( NumLib\\DomainDecomposition\\OGS5 FILES ${NumLib_DomainDecomposition_OGS5_Files} )


ADD_LIBRARY( NumLib STATIC 
	${NumLib_Files} 
	${NumLib_Core_Files} 
	${NumLib_TimeStepping_Files} 
	${NumLib_Coupling_Files} 
	${NumLib_Output_Files} 
	${NumLib_Discrete_Files} 
	${NumLib_DomainDecomposition_Files} 
	${NumLib_DomainDecomposition_OGS5_Files} 
)

INCLUDE_DIRECTORIES(
	.
	..
)

target_link_libraries (
	NumLib
	Base
        MathLib
        MeshLib
        FemLib
)
