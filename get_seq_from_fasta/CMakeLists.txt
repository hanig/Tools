cmake_minimum_required(VERSION 2.8.4)
project(get_seq_from_fasta)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/Modules/")

set(SOURCE_FILES main.c)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

add_executable(get_seq_from_fasta ${SOURCE_FILES})

find_package(Samtools REQUIRED)
include_directories(${Samtools_INCLUDE_DIRS})
set(LIBS ${LIBS} ${Samtools_LIBRARIES})

target_link_libraries(get_seq_from_fasta ${LIBS})