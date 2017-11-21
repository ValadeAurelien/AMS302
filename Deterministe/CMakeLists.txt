cmake_minimum_required(VERSION 3.5.2)

set(CMAKE_CXX_COMPILER g++)
project(deterministe)

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-O3") #"-O0 -g")

# Includes directory
include_directories(${PROJECT_SOURCE_DIR} /usr/local/include/Eigen)

# Sources directories
add_executable(solver_deter solver_deter.cpp legendre_polynomial.cpp)

