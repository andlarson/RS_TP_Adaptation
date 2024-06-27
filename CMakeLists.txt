cmake_minimum_required(VERSION 3.26 FATAL_ERROR)
project(RS_TP_Optimization LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(rs_tp_adaptation main.cxx)
target_include_directories(rs_tp_adaptation PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src/")
