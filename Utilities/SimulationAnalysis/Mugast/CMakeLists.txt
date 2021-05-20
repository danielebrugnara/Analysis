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

set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED OFF)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost 1.73 COMPONENTS serialization)

#if(Boost_FOUND)
include_directories(${Boost_INCLUDE_DIRS})
#add_executable(progname file1.cxx file2.cxx)
#target_link_libraries(progname ${Boost_LIBRARIES})
#endif()

#Adding the required root flags
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED)
include(${ROOT_USE_FILE})

include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_CXX_FLAGS})

include_directories($ENV{NPTOOL}/NPLib/include ${CMAKE_BINARY_DIR})
link_directories($ENV{NPTOOL}/NPLib/lib)

include_directories("${CMAKE_SOURCE_DIR}/../../../Lib/inc/")
#link_directories("${CMAKE_BINARY_DIR}/../../Build/Lib/bin/")
link_directories("${CMAKE_SOURCE_DIR}/../../../build/Lib/bin/")
link_directories("${CMAKE_BINARY_DIR}")


get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
	  message(STATUS "selected include dir='${dir}'")
endforeach()

get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY LINK_DIRECTORIES)
foreach(dir ${dirs})
	  message(STATUS "selected linking dir='${dir}'")
endforeach()

#static library for selector
ROOT_GENERATE_DICTIONARY(G__Selector Selector.h LINKDEF SelectorLinkDef.h)
add_library(Selector STATIC Selector.cxx G__Selector.cxx)
target_link_libraries(Selector ${ROOT_LIBRARIES} NPCore NPMugast NPMUST2 NPCATS NPInitialConditions NPReactionConditions NPInteractionCoordinates NPPhysics ReactionReconstruction Units EnergyLoss)

ROOT_GENERATE_DICTIONARY(G__SelectorData SelectorData.h LINKDEF SelectorDataLinkDef.h)
add_library(SelectorData STATIC SelectorData.cxx G__SelectorData.cxx)
target_link_libraries(SelectorData ${ROOT_LIBRARIES} MugastData Selector)

ROOT_GENERATE_DICTIONARY(G__RunSelector RunSelector.h LINKDEF RunSelectorLinkDef.h)
add_library(RunSelector STATIC RunSelector.cxx G__RunSelector.cxx)
target_link_libraries(RunSelector ${ROOT_LIBRARIES} ${Boost_LIBRARIES} Selector SelectorData)

#linking to main
add_executable(main main.cxx)
target_link_libraries(main RunSelector)
