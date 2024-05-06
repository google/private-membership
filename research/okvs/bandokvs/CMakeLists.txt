# main target
set(HEADER_FILES band.h band_okvs.h utils.h)
set(SRCS band.cc band_okvs.cc utils.cc)

add_library(bandokvs STATIC ${HEADER_FILES} ${SRCS})
target_link_libraries(bandokvs oc::libOTe)
#target_link_libraries(bandokvs libdivide)

if (APPLE)
    target_compile_options(bandokvs PUBLIC "-std=c++17")
else ()
    target_compile_options(bandokvs PUBLIC -g -std=c++17 -O3 -march=native)
endif ()
