include(${CMAKE_CURRENT_LIST_DIR}/preamble.cmake)

message(STATUS "OKVS_THIRDPARTY_DIR=${OKVS_THIRDPARTY_DIR}")


set(PUSHED_CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH})
set(CMAKE_PREFIX_PATH "${OKVS_THIRDPARTY_DIR};${CMAKE_PREFIX_PATH}")


#######################################
# sparsehash

#######################################
# libOTe

macro(FIND_LIBOTE)
    set(ARGS ${ARGN})

    #explicitly asked to fetch libOTe
    if (FETCH_LIBOTE)
        list(APPEND ARGS NO_DEFAULT_PATH PATHS ${OKVS_THIRDPARTY_DIR})
    endif ()

    find_package(libOTe ${ARGS})

    if (TARGET oc::libOTe)
        set(libOTe_FOUND ON)
    else ()
        set(libOTe_FOUND OFF)
    endif ()
endmacro()

if (FETCH_LIBOTE_AUTO)
    FIND_LIBOTE(QUIET)
    include(${CMAKE_CURRENT_LIST_DIR}/../thirdparty/getLibOTe.cmake)
endif ()

FIND_LIBOTE(REQUIRED)

macro(FIND_LIBDIVIDE)
    set(ARGS ${ARGN})

    #explicitly asked to fetch libdivide
    if(FETCH_LIBDIVIDE)
        list(APPEND ARGS NO_DEFAULT_PATH PATHS ${OKVS_THIRDPARTY_DIR})
    endif()

    find_path(LIBDIVIDE_INCLUDE_DIRS "libdivide.h" PATH_SUFFIXES "include" ${ARGS})
    if(EXISTS "${LIBDIVIDE_INCLUDE_DIRS}/libdivide.h")
        set(LIBDIVIDE_FOUND ON)
    else()
        set(LIBDIVIDE_FOUND OFF)
    endif()

endmacro()

#if(FETCH_LIBDIVIDE_AUTO)
#    FIND_LIBDIVIDE(QUIET)
#    include(${CMAKE_CURRENT_LIST_DIR}/../thirdparty/getLibDivide.cmake)
#endif()

#FIND_LIBDIVIDE(REQUIRED)

add_library(libdivide INTERFACE IMPORTED)

target_include_directories(libdivide INTERFACE
                $<BUILD_INTERFACE:${LIBDIVIDE_INCLUDE_DIRS}>
                $<INSTALL_INTERFACE:>)

# resort the previous prefix path
set(CMAKE_PREFIX_PATH ${PUSHED_CMAKE_PREFIX_PATH})
