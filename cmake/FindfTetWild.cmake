#[=======================================================================[.rst:
FindfTetWild
-------

Notes:
    Finds the fTetWild static library and associated include files so that it
        can be linked against.

Imported Targets:
    This module provides the following imported targets, if found:

    fTetWild::fTetWild
        The fTetWild static library.

Result Variables:
    This will define the following variables:
    
    fTetWild_FOUND
        True if the system has the fTetWild static library.
    fTetWild_INCLUDE_DIR
        Absolute path to include directory for fTetWild. 
    fTetWild_LIBRARY
        Absolute path to static library for fTetWild.

#]=======================================================================]

set(fTetWild_BUILD_DIR "/home/andlars/Downloads/fTetWild/build")

find_path(fTetWild_INCLUDE_DIR 
          NAMES FloatTetwild.h 
          PATHS "${fTetWild_BUILD_DIR}/include/floattetwild/" 
          NO_DEFAULT_PATH
         )

find_library(fTetWild_LIBRARY 
             NAMES libFloatTetwild.a 
             PATHS "${fTetWild_BUILD_DIR}" 
             NO_DEFAULT_PATH
            )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(fTetWild
                                  FOUND_VAR fTetWild_FOUND
                                  REQUIRED_VARS fTetWild_LIBRARY fTetWild_INCLUDE_DIR
                                 )

add_library(fTetWild::fTetWild UNKNOWN IMPORTED)
set_target_properties(fTetWild::fTetWild PROPERTIES
                      IMPORTED_LOCATION "${fTetWild_LIBRARY}"
                      INTERFACE_INCLUDE_DIRECTORIES "${fTetWild_INCLUDE_DIR}"
                     )

mark_as_advanced(fTetWild_INCLUDE_DIR fTetWild_LIBRARY fTetWild_FOUND)








