cmake_minimum_required(VERSION 2.8.4)
project(split_reads_by_length)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.c mystrio.c)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

add_executable(split_reads_by_length ${SOURCE_FILES})

find_package( ZLIB REQUIRED )
if ( ZLIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
    target_link_libraries( split_reads_by_length ${ZLIB_LIBRARIES} )
endif( ZLIB_FOUND )
