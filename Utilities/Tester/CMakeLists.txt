#CMakeLists.txt for analysis package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(tester)

#Where to find the libraries
#set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/Lib/bin")
#Where to find the executables
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}")

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

include_directories(${CMAKE_SOURCE_DIR}/../../Lib/inc ${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_CXX_FLAGS})

include_directories(${CMAKE_SOURCE_DIR}/../../Lib/inc)
include_directories(${CMAKE_SOURCE_DIR})
link_directories(${CMAKE_SOURCE_DIR})
link_directories(${CMAKE_SOURCE_DIR}/../../build/Lib/bin)


get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
      message(STATUS "selected include dir='${dir}'")
endforeach()

get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY LINK_DIRECTORIES)
foreach(dir ${dirs})
      message(STATUS "selected linking dir='${dir}'")
endforeach()

#link_directories(${CMAKE_SOURCE_DIR}/../../cmake-build-debug/Lib/bin)
add_executable(Tester Tester.cpp)
target_link_libraries(Tester ReactionReconstruction ReactionFragment Units ${ROOT_LIBRARIES})
