cmake_minimum_required(VERSION 3.8)
project(pyscan)

set(CMAKE_CXX_STANDARD 14)

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
        src/pyscan.cpp
        src/RectangleScan.cpp
        src/RectangleScan.hpp
        src/Statistics.hpp
        src/Utilities.hpp
        src/HalfplaneScan.cpp
        src/HalfplaneScan.hpp src/InOutPoint.cpp src/InOutPoint.hpp)

add_library(pyscan SHARED ${SOURCE_FILES})
