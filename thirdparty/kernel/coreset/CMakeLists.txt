cmake_minimum_required(VERSION 3.4.3)

project(coreset)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 14)


set(SOURCE_FILES
        appext.cpp
        generator.cpp
        myinc.cpp
        )


add_library(appext SHARED ${SOURCE_FILES})
target_link_libraries(appext ann)
