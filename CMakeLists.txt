cmake_minimum_required(VERSION 3.5)
project(pyscan)

set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_C_FLAGS_DEBUG "-DDEBUG")
find_package( Boost 1.36.0 COMPONENTS python REQUIRED)
find_package( PythonLibs 2.7 REQUIRED )
include_directories( ${PYTHON_INCLUDE_DIRS} )
include_directories(${Boost_INCLUDE_DIRS})
set(SOURCE_FILES
        src/BloomFilter.cpp
        src/BloomFilter.hpp
        src/EpsSamples.cpp
        src/EpsSamples.hpp
        src/FunctionApprox.cpp
        src/FunctionApprox.hpp
        src/Point.hpp
        src/RectangleScan.cpp
        src/RectangleScan.hpp
        src/Statistics.hpp
        src/Utilities.hpp
        src/HalfplaneScan.cpp
        src/HalfplaneScan.hpp
        src/DiskScan.cpp
        src/DiskScan.hpp
        src/Utilities.cpp src/DiskScan2.hpp src/DiskScan2.cpp)

set(PY_SOURCE_FILES
        src/pyscan.cpp src/DiskScan2.hpp src/DiskScan2.cpp)

set(test true)


if (test)
    set(TEST_FILES test/Utilities.hpp test/Utilities.cpp src/DiskScan2.hpp src/DiskScan2.cpp)

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

    # Now simply link against gtest or gtest_main as needed. Eg
    add_executable(DiskScan_test test/DiskScan_unittest.cpp src/HalfplaneScan.cpp src/DiskScan.cpp src/Utilities.cpp test/Utilities.cpp src/DiskScan2.cpp)
    target_link_libraries(DiskScan_test gtest_main)
    add_test(NAME DiskScan_test COMMAND DiskScan_test)

    add_executable(RectangleScan_test test/RectangleScan_unittest.cpp ${TEST_FILES} ${SOURCE_FILES})
    target_link_libraries(RectangleScan_test gtest_main)
    add_test(NAME RectangleScan_test COMMAND RectangleScan_test)

    add_executable(FunctionApprox_test test/FunctionApprox_unittest.cpp src/Utilities.cpp src/FunctionApprox.cpp test/Utilities.cpp)
    target_link_libraries(FunctionApprox_test gtest_main)
    add_test(NAME FunctionApprox_test COMMAND FunctionApprox_test)

    add_executable(Point_test test/Point_unittest.cpp)
    target_link_libraries(Point_test gtest_main)
    add_test(NAME Point_test COMMAND Point_test)

    add_executable(Halfplane_test test/Halfplane_unittest.cpp src/HalfplaneScan.cpp test/Utilities.cpp)
    target_link_libraries(Halfplane_test gtest_main)
    add_test(NAME Halfplane_test COMMAND Halfplane_test)
endif()


add_library(pyscan SHARED ${SOURCE_FILES} ${PY_SOURCE_FILES} src/DiskScan2.hpp src/DiskScan2.cpp)
target_link_libraries(pyscan ${Boost_LIBRARIES})
set_target_properties(pyscan PROPERTIES PREFIX "")
