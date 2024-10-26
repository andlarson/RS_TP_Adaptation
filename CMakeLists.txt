cmake_minimum_required(VERSION 3.30 FATAL_ERROR)
project(RS_TP_Optimization LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# -------------------------------------------------------------------------
#                               Dependencies
# -------------------------------------------------------------------------

# *** Geogram ***
# When installed, Geogram provides a find module. The find module still needs
#     a hint about where Geogram was installed, though.
# When I built and installed Geogram, the fresh version of GLFW3 was built
#     and installed. It lives at the installation location of Geogram. However,
#     the find module happens to look for system-installs of GLFW3 before looking
#     in the Geogram install directory. For this reason the system-install of
#     GLFW3 is used.
# list(APPEND CMAKE_MODULE_PATH "/Users/andrewlarson/Desktop/misc/programming/misc_installs/geogram_debug/lib/cmake/modules")
# set(GEOGRAM_INSTALL_PREFIX "/Users/andrewlarson/Desktop/misc/programming/misc_installs/geogram_debug/")
# find_package(Geogram REQUIRED MODULE)

# *** CGAL ***
find_package(CGAL REQUIRED COMPONENTS Qt6)

# -------------------------------------------------------------------------
#                               This Library     
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
                     )

# -------------------------------------------------------------------------
#                                   Tests 
# -------------------------------------------------------------------------

add_executable(test "${CMAKE_SOURCE_DIR}/test/test.cpp")

target_link_libraries(test
                      rs_tp_adaptation
                     )
