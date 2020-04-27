#CMakeLists.txt for analysis package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(simulationanalysis)

#Where to find the libraries
set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}")
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

include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_CXX_FLAGS})

include_directories($ENV{NPTOOL}/NPLib/include ${CMAKE_BINARY_DIR})
link_directories($ENV{NPTOOL}/NPLib/lib)

#static library for selector
ROOT_GENERATE_DICTIONARY(G__Selector Selector.h LINKDEF SelectorLinkDef.h)
add_library(Selector STATIC Selector.cxx G__Selector.cxx)
target_link_libraries(Selector ${ROOT_LIBRARIES} NPMugast NPMUST2 NPCATS NPInitialConditions NPReactionConditions NPInteractionCoordinates)

ROOT_GENERATE_DICTIONARY(G__RunSelector RunSelector.h LINKDEF RunSelectorLinkDef.h)
add_library(RunSelector STATIC RunSelector.cxx G__RunSelector.cxx)
target_link_libraries(RunSelector ${ROOT_LIBRARIES} Selector)

#linking to main
add_executable(main main.cxx)
target_link_libraries(main RunSelector)
