######################
### Find tools     ###
######################

# Find dot tool from graphviz
FIND_PROGRAM(DOT_TOOL_PATH dot DOC "Dot tool from graphviz")

# Find doxygen
FIND_PACKAGE(Doxygen)
#IF (DOXYGEN_FOUND)
#INCLUDE(scripts/cmake/cmake/DoxygenTargets.cmake)
#add_doxygen(scripts/cmake/cmake/DoxygenTargets.doxyfile.in)
#ENDIF(DOXYGEN_FOUND)

# Find gnu profiler gprof
FIND_PROGRAM(GPROF_PATH gprof DOC "GNU profiler gprof")

FIND_PACKAGE(cppcheck)


######################
### Find libraries ###
######################

FIND_PACKAGE(OpenMP)
IF(OPENMP_FOUND)
	SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
ENDIF()

