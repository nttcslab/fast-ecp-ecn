cmake_minimum_required(VERSION 3.8)
project(reliabilitylib)

set(SRCS ../beam_search/beam_search.cpp graph.cpp graphsimplify.cpp lwdp.cpp csrel.cpp)

add_library(reliabilitylib STATIC ${SRCS})

target_compile_features(reliabilitylib PUBLIC cxx_std_11)
#target_include_directories(reliabilitylib INTERFACE .)