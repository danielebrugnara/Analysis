#CMakeLists.txt for analysis package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(effanalysis)

#Where to find the libraries
set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}")
#Where to find the executables
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}")

#For a debug build compile with cmake .. -DCMAKE_BUILD_TYPE=Debug
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
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

include_directories("${CMAKE_BINARY_DIR}/../../Lib/inc/")
include_directories("${CMAKE_BINARY_DIR}")
link_directories("${CMAKE_BINARY_DIR}/../../Build/Lib/bin/")
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
target_link_libraries(Selector ${ROOT_LIBRARIES})

ROOT_GENERATE_DICTIONARY(G__RunSelector RunSelector.h LINKDEF RunSelectorLinkDef.h)
add_library(RunSelector STATIC RunSelector.cxx G__RunSelector.cxx)
target_link_libraries(RunSelector ${ROOT_LIBRARIES} Selector)

add_library(Units STATIC Units.h)
set_target_properties(Units PROPERTIES LINKER_LANGUAGE CXX)

add_library(DiaGraph STATIC DiaGraph.h)
set_target_properties(DiaGraph PROPERTIES LINKER_LANGUAGE CXX)

add_library(LevelScheme STATIC LevelScheme.h LevelScheme.cxx)
target_link_libraries(LevelScheme DiaGraph Units)

ROOT_GENERATE_DICTIONARY(G__SpectrumAnalyzer SpectrumAnalyzer.h LINKDEF SpectrumAnalyzerLinkDef.h)
add_library(SpectrumAnalyzer STATIC SpectrumAnalyzer.cxx G__SpectrumAnalyzer.cxx)
target_link_libraries(SpectrumAnalyzer ${ROOT_LIBRARIES} DiaGraph LevelScheme)

#linking to main
add_executable(effanalysis main.cxx)
target_link_libraries(effanalysis RunSelector SpectrumAnalyzer)
