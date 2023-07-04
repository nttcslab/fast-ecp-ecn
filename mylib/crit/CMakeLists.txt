cmake_minimum_required(VERSION 3.8)
project(reliability-critlib)

set(SRCS ../../beam_search/beam_search.cpp graph.cpp attr.cpp)

add_library(reliability-critlib STATIC ${SRCS})

find_library(ADEPTLIB NAMES adept)
target_link_libraries(reliability-critlib ${ADEPTLIB})
target_compile_features(reliability-critlib PUBLIC cxx_std_11)
#target_include_directories(reliability-critlib INTERFACE .)