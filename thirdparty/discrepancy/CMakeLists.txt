cmake_minimum_required(VERSION 3.4.3)
project(discrepancy)



include_directories(include/)

set(SOURCE_FILES
        src/dataset.c
        src/find_rect.c
        src/function.c
        src/interval.c
        src/mtint.c
        src/naive_scan.c
        src/naive_scan_grid.c
        src/scan_grid.c
        src/scan_stat.c
        src/sort.c
        src/stat_rect.c)
#set(HEADER ${CMAKE_SOURCE_DIR}/include/*.h)


add_library(discrepancy SHARED ${SOURCE_FILES})
#target_link_libraries(libpyscan PRIVATE ${THIRDPARTY_LIBRARIES})
