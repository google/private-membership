cmake_minimum_required(VERSION 3.15)

set(CMAKE_CXX_COMPILER clang++-14)

project(compression VERSION 0.1.0)

set(COMPRESSION_BUILD ON)

# sets some variables
include(cmake/preamble.cmake)

# sets build option
include(cmake/buildOptions.cmake)

# find all the dependancies.
include(cmake/findDependancies.cmake)

add_subdirectory(bandokvs)
add_subdirectory(frontend)

target_link_libraries(frontend bandokvs)

# setup the install
include(cmake/install.cmake)