cmake_minimum_required(VERSION 3.30 FATAL_ERROR)
project(RS_TP_Optimization LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# -------------------------------------------------------------------------
#                               Dependencies
# -------------------------------------------------------------------------

# When installed, Geogram provides a find module.
list(APPEND CMAKE_MODULE_PATH "/Users/andrewlarson/Desktop/misc/programming/misc_installs/geogram_debug/lib/cmake/modules")
find_package(Geogram REQUIRED MODULE)

# When I built Geogram, I set it up to also build GLFW. When installed,
#     Geogram also contains GLFW. GLFW provides a CMake package configuration 
#     file.
# set(glfw3_ROOT "/Users/andrewlarson/Desktop/misc/programming/misc_installs/geogram_debug/lib/cmake/glfw3/")
# find_package(glfw3 REQUIRED CONFIG)


# -------------------------------------------------------------------------
#                               This Library     
# -------------------------------------------------------------------------

add_library(rs_tp_adaptation SHARED
            "${CMAKE_SOURCE_DIR}/src/main.cpp"
           )

target_link_libraries(rs_tp_adaptation 
                      Geogram::geogram
                      Geogram::geogram_gfx
                     )

# -------------------------------------------------------------------------
#                                   Tests 
# -------------------------------------------------------------------------
