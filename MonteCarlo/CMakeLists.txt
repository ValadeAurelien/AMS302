cmake_minimum_required(VERSION 3.5.2)

set(CMAKE_CXX_COMPILER g++)
project(MonteCarlo)

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

set (CMAKE_CXX_STANDARD 11)
set(CXX_FLAGS "-Wall -O3")

add_executable(TP1 TP1.cpp)
add_executable(TP1_sigmavar TP1_sigmavar.cpp)
add_executable(TP1_diff TP1_diff.cpp)
