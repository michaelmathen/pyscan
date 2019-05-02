cmake_minimum_required(VERSION 3.4.3)

project(pyscan)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 17)

#Need this for GSL
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_FLAGS_DEBUG "-fPIC -Wall -Wextra -g -O1 -fno-omit-frame-pointer")
set(CMAKE_CXX_FLAGS_RELEASE "-fPIC -w -O2 -march=native -DNDEBUG")
set(CMAKE_C_FLAGS_DEBUG "-fPIC -Wall -Wextra -g -O1 -fno-omit-frame-pointer")
set(CMAKE_C_FLAGS_RELEASE "-fPIC -w -O2 -march=native -DNDEBUG")

add_subdirectory(thirdparty/pybind11)


find_package(Boost REQUIRED)
find_package(CGAL REQUIRED)
find_package(GSL REQUIRED)

#include( ${CGAL_USE_FILE} )
#set(CGAL_DONT_OVERRIDE_CMAKE_FLAGS TRUE CACHE BOOL "Force CGAL to maintain CMAKE flags")
#include(${CGAL_USE_FILE})


execute_process(COMMAND ${CMAKE_COMMAND} .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/thirdparty/kernel/ANN/ )
if(result)
    message(FATAL_ERROR "CMake step for ANN: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/thirdparty/kernel/ANN/ )
if(result)
    message(FATAL_ERROR "build failed for ANN: ${result}")
endif()

execute_process(COMMAND ${CMAKE_COMMAND} .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/thirdparty/kernel/coreset/ )
if(result)
    message(FATAL_ERROR "CMake step for coreset: ${result}")
endif()


execute_process(COMMAND ${CMAKE_COMMAND} --build .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/thirdparty/kernel/coreset )
if(result)
    message(FATAL_ERROR "build failed for coreset: ${result}")
endif()

add_subdirectory(thirdparty/discrepancy)
set(THIRDPARTY_LIBRARIES
        discrepancy
        #${CMAKE_BINARY_DIR}/thirdparty/discrepancy/libdiscrepancy.so
        ${CMAKE_SOURCE_DIR}/thirdparty/kernel/coreset/libappext.so
        ${CMAKE_SOURCE_DIR}/thirdparty/kernel/ANN/libann.so
        CGAL::CGAL
        GSL::gsl
        GSL::gslcblas)

include_directories(
        ${Boost_INCLUDE_DIRS}
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/thirdparty/discrepancy/include
        ${CMAKE_SOURCE_DIR}/thirdparty/kernel/coreset/include
        ${CMAKE_SOURCE_DIR}/thirdparty/discrepancy/include
)

set(SOURCE_FILES
        src/BloomFilter.cpp
        src/FunctionApprox.cpp
        src/RectangleScan.cpp
        src/HalfSpaceScan.cpp
        src/DiskScan.cpp
        src/Point.cpp
        src/Utilities.cpp
        src/Test_Utilities.cpp
        src/TrajectoryScan.cpp
        src/TrajectoryCoreSet.cpp
        src/ConvexHull.cpp
        src/Segment.cpp
        src/HamTree.cpp
        src/IntervalScan.cpp
        src/RegionCoreSet.cpp
        src/JeffCodes.cpp
        src/AnnuliScanning.cpp
        #src/kernel.cpp
        )

#set(SOURCE_JEFF_FILES
#        )

set(HEADER include/BloomFilter.hpp
        include/DiskScan.hpp
        include/FunctionApprox.hpp
        include/HalfSpaceScan.hpp
        include/Point.hpp
        include/RectangleScan.hpp
        include/RegionCoreSet.hpp
        include/Segment.hpp
        include/SparseGrid.hpp
        include/Statistics.hpp
        include/Test_Utilities.hpp
        include/ConvexHull.hpp
        include/TrajectoryScan.hpp
        include/Sampling.hpp
        include/HamTree.hpp
        inlcude/AnnuliScanning.hpp
        include/IntervalScan.hpp
        include/Utilities.hpp)

set(PY_SOURCE_FILES
        src/pyscan.cpp)


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

SET(TEST_LINK_LIBRARIES gtest_main ${THIRDPARTY_LIBRARIES})

#add_library(discrepancy_obj OBJECT ${SOURCE_JEFF_FILES})
add_library(src_obj OBJECT ${SOURCE_FILES})
########################################
#Compile various tests##################
########################################

set(TEST_NAMES
        DiskScan_unittest
        FunctionApprox_unittest
        Halfplane_unittest
        Testing
        Point_unittest
        RectangleScan_unittest
        trajectory_unittest)


foreach(TEST ${TEST_NAMES})

    add_executable(${TEST} ${CMAKE_SOURCE_DIR}/test/${TEST}.cpp $<TARGET_OBJECTS:src_obj>)
    target_link_libraries(${TEST} ${TEST_LINK_LIBRARIES} )
    add_test(NAME ${TEST} COMMAND ${TEST})

endforeach()


#######################################################
#Compile the python wrapper############################
#######################################################
configure_file(pyscan.py ${CMAKE_CURRENT_BINARY_DIR}/pyscan.py COPYONLY)
pybind11_add_module(libpyscan SHARED $<TARGET_OBJECTS:src_obj> ${PY_SOURCE_FILES})
#add_library(pyscan SHARED $<TARGET_OBJECTS:src_obj> ${PY_SOURCE_FILES})
target_link_libraries(libpyscan PRIVATE ${THIRDPARTY_LIBRARIES})
