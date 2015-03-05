cmake_minimum_required(VERSION 2.8.4)
project(generate_count_table_bamfile)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/Modules/")

set(SOURCE_FILES main.cpp sequences.c dataio.c statistics.c hashtable.c information.c interval_tree.cpp readio.c)
add_executable(generate_count_table_bamfile ${SOURCE_FILES})


find_package(Samtools REQUIRED)
include_directories(${Samtools_INCLUDE_DIRS})
set(LIBS ${LIBS} ${Samtools_LIBRARIES})

target_link_libraries(generate_count_table_bamfile ${LIBS})