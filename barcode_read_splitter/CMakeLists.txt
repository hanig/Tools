cmake_minimum_required(VERSION 2.8.4)
project(barcode_read_splitter)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.c mystrio.c hashtable.c)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

add_executable(barcode_read_splitter ${SOURCE_FILES})

find_package( ZLIB REQUIRED )
if ( ZLIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
    target_link_libraries( barcode_read_splitter ${ZLIB_LIBRARIES} )
endif( ZLIB_FOUND )
