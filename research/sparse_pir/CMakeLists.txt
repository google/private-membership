# CMakeList.txt : Top-level CMake project file, do global configuration
# and include sub-projects here.
#
cmake_minimum_required (VERSION 3.8)

project ("sparse_pir")

# Include sub-projects.

add_subdirectory ("Onion-PIR")
add_subdirectory ("sparse_pir")
