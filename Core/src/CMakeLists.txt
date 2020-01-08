include_directories(${analysis_SOURCE_DIR}/Core/inc)
include_directories(${analysis_SOURCE_DIR}/Lib/inc)
include_directories($ENV{NPTOOL}/NPLib/include)


link_directories(${CMAKE_BINARY_DIR}/Lib/bin)
link_directories($ENV{NPTOOL}/NPLib/lib)

add_executable(analysis Main.cxx)
target_link_libraries(analysis Analysis)