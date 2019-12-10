#CMakeLists.txt for analysis package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(analysis)

### where to find the libraries
set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/Lib/bin")
### where to find the executables
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED)
include(${ROOT_USE_FILE})

include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_CXX_FLAGS})


add_subdirectory(Lib)
add_subdirectory(Core)

install(DIRECTORY ${CMAKE_SOURCE_DIR}/Configs DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

