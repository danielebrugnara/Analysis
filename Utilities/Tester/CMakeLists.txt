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
#set(CMAKE_CXX_FLAGS_DEBUG "-g")
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

#include_directories(${CMAKE_SOURCE_DIR}/../../Lib/inc ${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_CXX_FLAGS})

add_library(Units SHARED Units.cxx)
set_target_properties(Units PROPERTIES PUBLIC_HEADER ${analysis_SOURCE_DIR}/Lib/inc)

ROOT_GENERATE_DICTIONARY(G__ReactionFragment ReactionFragment.h LINKDEF ReactionFragmentLinkDef.h)
add_library(ReactionFragment SHARED ReactionFragment.cxx G__ReactionFragment.cxx)
target_link_libraries(ReactionFragment ${ROOT_LIBRARIES})

ROOT_GENERATE_DICTIONARY(G__ReactionReconstruction ReactionReconstruction.h LINKDEF ReactionReconstructionLinkDef.h)
add_library(ReactionReconstruction SHARED ReactionReconstruction.cxx G__ReactionReconstruction.cxx)
target_link_libraries(ReactionReconstruction ${ROOT_LIBRARIES} ReactionFragment)

#include_directories(${CMAKE_SOURCE_DIR}/../../Lib/inc)
include_directories(${CMAKE_SOURCE_DIR})
link_directories(${CMAKE_SOURCE_DIR})
#link_directories(${CMAKE_SOURCE_DIR}/../../cmake-build-debug/Lib/bin)
add_executable(Tester Tester.cpp)
target_link_libraries(Tester ReactionReconstruction ReactionFragment Units ${ROOT_LIBRARIES})
