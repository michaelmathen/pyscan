cmake_minimum_required(VERSION 3.4.3)

project(coreset)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -g")
set(CMAKE_CXX_FLAGS_RELEASE "-fPIC -w -g")


set(SOURCE_FILES
        appext.cpp
        generator.cpp
        myinc.cpp
        )


set(HEADER include/appext.h
        include/generator.h
        include/myinc.h)

SET(ANN_LIB ${CMAKE_BINARY_DIR}/../ANN/libann.so)

include_directories( ${CMAKE_SOURCE_DIR}/include ${CMAKE_BINARY_DIR}/../ANN/include/)

add_library(appext SHARED ${SOURCE_FILES})
target_link_libraries(appext ${ANN_LIB})
