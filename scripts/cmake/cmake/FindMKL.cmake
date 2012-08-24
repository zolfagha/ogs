# Find Intel Math Karnel Library

#if (MKL_INCLUDES AND MKL_LIBRARIES)
#  # Already in cache, be silent
#  set (MKL_FIND_QUIETLY TRUE)
#endif (MKL_INCLUDES AND MKL_LIBRARIES)

FIND_PATH(MKL_INCLUDES NAMES mkl.h
    HINTS ENV MKL_DIR 
)

# Tell if the unix system is on 64-bit base
if (${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "x86_64")
#if(CMAKE_SIZEOF_VOID_P MATCHES "8")
    MESSAGE (STATUS  "Look for 64-bit MKL" )
    find_library(MKL_LIBRARIES
      mkl_core
      PATHS
      $ENV{MKLLIB}
            /opt/intel/mkl/*/lib/em64t
            /opt/intel/Compiler/*/*/mkl/lib/em64t
      ${LIB_INSTALL_DIR}
    )

    if(MKL_LIBRARIES)
        set(MKL_LIBRARIES ${MKL_LIBRARIES} mkl_intel_lp64 mkl_solver_lp64)
        if(CMAKE_C_COMPILER EQUAL "icc")
            set(MKL_LIBRARIES ${MKL_LIBRARIES} mkl_intel_thread iomp5 pthread)
        else()
            set(MKL_LIBRARIES ${MKL_LIBRARIES} mkl_gnu_thread)
        endif()
    endif()

else()
    MESSAGE (STATUS  "Look for 32-bit MKL" )
    find_library(MKL_LIBRARIES
      mkl_core
      PATHS
      $ENV{MKLLIB}
      /opt/intel/mkl/*/lib/32
      /opt/intel/Compiler/*/*/mkl/lib/32
      ${LIB_INSTALL_DIR}
    )
    
    if(MKL_LIBRARIES)
        set(MKL_LIBRARIES ${MKL_LIBRARIES} mkl_intel mkl_solver)
        if(CMAKE_C_COMPILER EQUAL "icc")
            set(MKL_LIBRARIES ${MKL_LIBRARIES} mkl_intel_thread iomp5 pthread)
        else()
            set(MKL_LIBRARIES ${MKL_LIBRARIES} mkl_gnu_thread)
        endif()
    endif()

endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MKL DEFAULT_MSG MKL_INCLUDES MKL_LIBRARIES)
mark_as_advanced(MKL_INCLUDES MKL_LIBRARIES)

