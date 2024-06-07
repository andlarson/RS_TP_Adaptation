cmake_minimum_required(VERSION 3.26 FATAL_ERROR)
project(RS_TP_Optimization LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Search path for the find modules.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")

# Use the find module!
find_package(fTetWild MODULE REQUIRED)

add_executable(rs_tp_adapation main.cxx)
target_link_libraries(rs_tp_adaptation fTetWild::fTetWild)
