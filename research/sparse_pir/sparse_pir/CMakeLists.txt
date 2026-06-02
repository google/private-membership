# CMakeList.txt : CMake project for sparse_pir, include source and define
# project specific logic here.
#
cmake_minimum_required(VERSION 3.8)

# Add source to this project's executable.
set(HEADER_FILES "band_row_vector.h" "bin_allocator.h" "band_row_vector_hasher.h" "int_utils.h" "prime_field.h" "linear_system.h" "linear_system_row.h" "main.cpp" "sparse_pir.h" "plaintext_value.h")
set(SRC_FILES "main.cpp")
add_executable(sparse_pir ${HEADER_FILES} ${SRC_FILES})

if (CMAKE_VERSION VERSION_GREATER 3.12)
    set_property(TARGET sparse_pir PROPERTY CXX_STANDARD 20)
endif ()

find_package(SEAL 3.5 REQUIRED)
target_link_libraries(sparse_pir SEAL::seal)
find_package(libOTe REQUIRED)
target_link_libraries(sparse_pir oc::libOTe)
target_compile_options(sparse_pir PUBLIC -msse4.2 -mavx -mavx2 -mfma -mavx512f -mavx512vl -mavx512dq)
