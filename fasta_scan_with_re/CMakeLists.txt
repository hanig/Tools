cmake_minimum_required(VERSION 2.8.4)
project(fasta_scan_with_re)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.c dataio.c sequences.c)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

add_executable(fasta_scan_with_re ${SOURCE_FILES})
SET_TARGET_PROPERTIES(fasta_scan_with_re
    PROPERTIES
    COMPILE_FLAGS "-DUSE_PCRE"
    LINK_FLAGS "-L/usr/local/lib -lpcre"
)