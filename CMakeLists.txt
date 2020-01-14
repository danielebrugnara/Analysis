#CMakeLists.txt for analysis package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(analysis)

#Where to find the libraries
set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/Lib/bin")
#Where to find the executables
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")

#For a debug build compile with cmake .. -DCMAKE_BUILD_TYPE=Debug
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

#Verbose option for debugging cmake .. -DVERBOSE_DEBUG=ON
option(VERBOSE_DEBUG "Add extra couts" OFF)

if(VERBOSE_DEBUG)
    add_definitions(-DVERBOSE_DEBUG)
endif()

#Adding the required root flags
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED)
include(${ROOT_USE_FILE})

include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_CXX_FLAGS})

#Dynamic libraries
add_subdirectory(Lib)

#Main executable source
add_subdirectory(Core)

install(DIRECTORY ${CMAKE_SOURCE_DIR}/Configs DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(DIRECTORY ${CMAKE_SOURCE_DIR}/Out DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

