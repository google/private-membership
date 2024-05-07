# these are just pass through config file for the ones that are placed in the build directory.


include("${CMAKE_CURRENT_LIST_DIR}/preamble.cmake")

if (NOT EXISTS "${OKVS_BUILD_DIR}")
    message(FATAL_ERROR "failed to find the okvs build directory. Looked at OKVS_BUILD_DIR: ${OKVS_BUILD_DIR}\n Please set it manually.")
endif ()

include("${OKVS_BUILD_DIR}/okvsConfig.cmake")