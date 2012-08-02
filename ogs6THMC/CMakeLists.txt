

# Scans the current directory and returns a list of subdirectories.
# Third parameter is 1 if you want relative paths returned.
# Usage: list_subdirectories(the_list_is_returned_here C:/cwd 1)
macro(list_subdirectories retval curdir return_relative)
#   MESSAGE (STATUS "curdir: ${curdir}")
#   file(GLOB sub-dir ${curdir}/*)
    file(GLOB sub-dir RELATIVE ${curdir} ${curdir}/*)
#   list(LENGTH sub-dir n_dir)
#   MESSAGE (STATUS "sub-dir (${n_dir}): ${sub-dir}")
    set(list_of_dirs "")
    foreach(dir ${sub-dir})
        if(IS_DIRECTORY ${curdir}/${dir})
            if (${return_relative})
                set(list_of_dirs ${list_of_dirs} ${dir})
            else()
                set(list_of_dirs ${list_of_dirs} ${curdir}/${dir})
            endif()
        endif()
    endforeach()
    set(${retval} ${list_of_dirs})
endmacro()


#look for sub directories in ProcessLib
#MESSAGE (STATUS "Look for subdirectory")
list_subdirectories(List_subdir ${CMAKE_CURRENT_SOURCE_DIR}/modules 1)
LIST(LENGTH List_subdir n_dir)
MESSAGE (STATUS "${n_dir} processses are found.")

#add soruce files
FOREACH(dir ${List_subdir})
#    MESSAGE (STATUS "sub dir found: ${dir}")
    IF(EXISTS modules/${dir}/CMakeLists.txt)
        INCLUDE(modules/${dir}/CMakeLists.txt)
    ENDIF()
    FILE(GLOB TempProcessLib_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} modules/${dir}/*.h modules/${dir}/*.cpp)
    SOURCE_GROUP( ogs6THMC\\modules\\${dir} FILES ${TempProcessLib_Files} )
    LIST(APPEND ProcessLib_Files ${TempProcessLib_Files}) 
ENDFOREACH()

FILE(GLOB TempProcessLib_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} modules/*.h modules/*.cpp modules/*.txt)
SOURCE_GROUP( ogs6THMC\\modules FILES ${TempProcessLib_Files} )
LIST(APPEND ProcessLib_Files ${TempProcessLib_Files}) 

FILE(GLOB TempProcessLib_Files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h *.cpp)
SOURCE_GROUP( ogs6THMC FILES ${TempProcessLib_Files} )
LIST(APPEND ProcessLib_Files ${TempProcessLib_Files}) 

ADD_LIBRARY( ogs6THMC STATIC 
    ${ProcessLib_Files}
)

SET_TARGET_PROPERTIES(ogs6THMC PROPERTIES LINKER_LANGUAGE CXX)

INCLUDE_DIRECTORIES(
        .
        ../ogsCore
        ../ogsGeometry
        ../ogsMesh
        ../ogsNumerics
        ../ogsPhysics
        ../ogsFileIO
        ../ogsProcess
        ../extra/Eigen
        ../extra/tinyxml2
)

target_link_libraries( ogs6THMC
    ogsCore
    ogsGeometry
    ogsMesh
    ogsNumerics
    ogsFileIO
    ogsProcess
    tinyxml2
)

SET_TARGET_PROPERTIES(ogs6THMC PROPERTIES LINKER_LANGUAGE CXX)

