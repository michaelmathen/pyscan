cmake_minimum_required(VERSION 3.4.3)

project(pyscan)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -g")
set(CMAKE_CXX_FLAGS_RELEASE "-fPIC -w -g -fsanitize=address")




find_package( PythonLibs 3 REQUIRED )
include_directories( ${PYTHON_INCLUDE_DIRS} )

set(version ${PYTHONLIBS_VERSION_STRING})
STRING( REGEX REPLACE "([0-9.]+).[0-9]+(rc[0-9]+)?" "\\1" version ${version} )
STRING( REGEX REPLACE "[^0-9]" "" version ${version} )
find_package(Boost COMPONENTS "python${version}" REQUIRED)


execute_process(COMMAND ${CMAKE_COMMAND} .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/thirdparty/libgdiam-1.0.3 )
if(result)
    message(FATAL_ERROR "CMake step for libgdiam failed: ${result}")
endif()


execute_process(COMMAND ${CMAKE_COMMAND} --build .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/thirdparty/libgdiam-1.0.3 )
if(result)
    message(FATAL_ERROR "Build step for libgdiam failed: ${result}")
endif()

set(THIRDPARTY_LIBRARIES ${CMAKE_SOURCE_DIR}/thirdparty/libgdiam-1.0.3/libgdiam.so)

include_directories(
        ${Boost_INCLUDE_DIRS}
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/thirdparty/libgdiam-1.0.3)

set(SOURCE_FILES
        src/BloomFilter.cpp
        src/FunctionApprox.cpp
        src/RectangleScan.cpp
        src/HalfSpaceScan.cpp
        src/DiskScan.cpp
        src/Point.cpp
        src/Utilities.cpp
        src/Test_Utilities.cpp
        src/DiskScan2.cpp
        src/TrajectoryScan.cpp
        #src/kernel.cpp
        )

set(HEADER include/BloomFilter.hpp
        include/DiskScan.hpp
        include/DiskScan2.hpp
        include/FunctionApprox.hpp
        include/HalfSpaceScan.hpp
        include/Point.hpp
        include/RectangleScan.hpp
        include/SparseGrid.hpp
        include/Statistics.hpp
        include/Test_Utilities.hpp
        include/TrajectoryScan.hpp
        include/Utilities.hpp)

set(PY_SOURCE_FILES
        src/pyscan.cpp src/Point.cpp)

set(TEST_FILES test/Test_Utilities.cpp ${SOURCE_FILES})

#################################################################
# Setup google test files.#######################################
#################################################################

configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
if(result)
    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
if(result)
    message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_BINARY_DIR}/googletest-src
        ${CMAKE_BINARY_DIR}/googletest-build)

# The gtest/gtest_main targets carry header search path
# dependencies automatically when using CMake 2.8.11 or
# later. Otherwise we have to add them here ourselves.
if (CMAKE_VERSION VERSION_LESS 2.8.11)
  include_directories("${gtest_SOURCE_DIR}/include")
endif()

SET(TEST_LINK_LIBRARIES gtest_main)

########################################
#Compile various tests##################
########################################

set(TEST_NAMES
        DiskScan_unittest
        FunctionApprox_unittest
        Halfplane_unittest
        Point_unittest
        RectangleScan_unittest
        trajectory_unittest)


foreach(TEST ${TEST_NAMES})

    add_executable(${TEST} ${CMAKE_SOURCE_DIR}/test/${TEST}.cpp ${SOURCE_FILES})
    target_link_libraries(${TEST} ${TEST_LINK_LIBRARIES})
    add_test(NAME ${TEST} COMMAND ${TEST})

endforeach()






#######################################################
#Compile the python wrapper############################
#######################################################
configure_file(pyscan.py ${CMAKE_CURRENT_BINARY_DIR}/pyscan.py COPYONLY)
add_library(pyscan SHARED ${SOURCE_FILES} ${PY_SOURCE_FILES})
target_link_libraries(pyscan ${Boost_LIBRARIES} ${THIRDPARTY_LIBRARIES})
