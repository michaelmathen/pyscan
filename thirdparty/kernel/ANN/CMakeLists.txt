cmake_minimum_required(VERSION 3.4.3)

project(ANN)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -g")
set(CMAKE_CXX_FLAGS_RELEASE "-fPIC -w -g  -ansi")


set(SOURCE_FILES
        src/bd_pr_search.cc
        src/ANN.cc
        src/bd_tree.cc
        src/bd_search.cc
        src/brute.cc
        src/kd_pr_search.cc
        src/kd_search.cc
        src/kd_split.cc
        src/kd_tree.cc
        src/kd_util.cc
        src/perf.cc
        )


add_library(ann SHARED ${SOURCE_FILES})
#target_link_libraries(pyscan)
