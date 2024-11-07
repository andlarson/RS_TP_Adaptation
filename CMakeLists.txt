cmake_minimum_required(VERSION 3.26 FATAL_ERROR)
project(RS_TP_Optimization LANGUAGES CXX)

# -------------------------------------------------------------------------
#                                Preamble 
# -------------------------------------------------------------------------

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include("${CMAKE_SOURCE_DIR}/cmake/utils.cmake")
GCC8_filesystem_stl_fix()

# -------------------------------------------------------------------------
#                               Dependencies
# -------------------------------------------------------------------------

# Danger: Config mode must be forced because Eigen provides a Find Module that
#     doesn't provide the Eigen component. It's unclear why the Find Module
#     exists when a Configuration file is provided.
# This is a very old version of Eigen. I will probably upgrade to a newer
#     version.
find_package(Eigen3 3.3 REQUIRED CONFIG)

# Danger: On some systems, version 6.0.1 of CGAL and its dependencies aren't
#     included in the system's package manager. However, it's very desirable
#     to have version 6.0.1 of CGAL because it adds significant visualization
#     functionality.
#
# Precisely, the required dependencies of CGAL 6.0.1 are:
#     1) A compiler that supports C++17.
#     2) CMake version 3.22 or newer.
#     3) Boost version 1.72 or newer.
#     4) ((GMP 5.0.1 or newer AND MPFR 3.0.0 or newer) OR (Boost version 1.72 or
#            newer)) for multiprecision arithmetic.
# The optional dependencies of CGAL 6.0.1 that this repository uses include:
#     1) Qt version 6.4 or newer.
#
# Most systems have:
#     1) A compiler that supports C++17.
# Of the dependenices, most system's package manager has:
#     1) CMake version 3.22 or newer.
#     2) GMP 5.0.1 or newer and MPFR 3.0.0 or newer.
# Of the dependencies, the system's package manager sometimes doesn't have:
#     1) Boost version 1.72 or newer. 
#     2) Qt version 6.4 or newer. 
# Therefore, I've often: 
#     1) Installed Boost 1.72 or newer from source.
#            CGAL only needs the header-only portion of the library.
#     2) Installed Qt version 6.4 or newer using the Qt Installer GUI program. 
#            This could be built from source in the future. 
#
# Danger: For some reason, when I installed and used Qt version 6.8.0, some
#     obscure error is generated when compiling a Qt? object file. To fix this,
#     I downgraded to Qt version 6.7.3.
find_package(CGAL 6.0.1 REQUIRED COMPONENTS Qt6)

# -------------------------------------------------------------------------
#                                 Library     
# -------------------------------------------------------------------------

add_library(rs_tp_adaptation SHARED
            "${CMAKE_SOURCE_DIR}/src/estimation.cpp"
            "${CMAKE_SOURCE_DIR}/src/mesh_utilities.cpp"
           )

target_include_directories(rs_tp_adaptation 
                           PUBLIC "include/"
                           PRIVATE "${CMAKE_SOURCE_DIR}/src/include")

target_link_libraries(rs_tp_adaptation 
                      CGAL::CGAL
                      CGAL::CGAL_Basic_viewer
                      Eigen3::Eigen
                     )

# -------------------------------------------------------------------------
#                                   Tests 
# -------------------------------------------------------------------------

add_executable(test "${CMAKE_SOURCE_DIR}/test/test.cpp")

target_compile_definitions(test PUBLIC "TEST_DIR=\"${CMAKE_SOURCE_DIR}/test/\"")

target_link_libraries(test
                      rs_tp_adaptation
                     )
