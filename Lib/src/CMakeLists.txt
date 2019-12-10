#include_directories(${analysis_SOURCE_DIR}/Lib/inc $ENV{NPTOOL}/NPLib/include ${ROOT_INCLUDE_DIRS})
include_directories(${analysis_SOURCE_DIR}/Lib/inc $ENV{NPTOOL}/NPLib/include)
#add_definitions(${ROOT_CXX_FLAGS})  

#add_library(mugast SHARED IMPORTED)
#set_property(TARGET mugast PROPERTY IMPORTED_LOCATION "$ENV{NPTOOL}/NPLib/libs/libNPMugast.so")
#
#add_library(must2 SHARED IMPORTED)
#set_property(TARGET mugast PROPERTY IMPORTED_LOCATION "$ENV{NPTOOL}/NPLib/libs/libNPCATS.so")
#
#add_library(cats SHARED IMPORTED)
#set_property(TARGET mugast PROPERTY IMPORTED_LOCATION "$ENV{NPTOOL}/NPLib/libs/libNPMUST2.so")

link_directories($ENV{NPTOOL}/NPLib/lib)


ROOT_GENERATE_DICTIONARY(G__Interpolation Interpolation.h LINKDEF ${analysis_SOURCE_DIR}/Lib/inc/InterpolationLinkDef.h)
add_library(Interpolation SHARED Interpolation.cxx G__Interpolation.cxx)

ROOT_GENERATE_DICTIONARY(G__Selector Selector.h LINKDEF ${analysis_SOURCE_DIR}/Lib/inc/SelectorLinkDef.h)
add_library(Selector SHARED Selector.cxx G__Selector.cxx)


ROOT_GENERATE_DICTIONARY(G__Analysis Analysis.h LINKDEF ${analysis_SOURCE_DIR}/Lib/inc/AnalysisLinkDef.h)
add_library(Analysis SHARED Analysis.cxx G__Analysis.cxx)

#ROOT_GENERATE_DICTIONARY(G__Selector Selector.h LINKDEF SelectorLinkDef.h)
#add_library(Selector SHARED Selector.cxx G__Selector.cxx)

#ROOT_GENERATE_DICTIONARY(G__RunAnalysis RunAnalysis.h LINKDEF RunAnalysisLinkDef.h)
#add_library(RunAnalysis SHARED RunAnalysis.cxx G__RunAnalysis.cxx)

#ROOT_GENERATE_DICTIONARY(G__Analysis Analysis.h LINKDEF AnalysisLinkDef.h)
#add_library(Analysis SHARED Analysis.cxx G__Analysis.cxx)

set_target_properties(Selector PROPERTIES PUBLIC_HEADER ${analysis_SOURCE_DIR}/Lib/inc)
set_target_properties(Interpolation PROPERTIES PUBLIC_HEADER ${analysis_SOURCE_DIR}/Lib/inc)
set_target_properties(Analysis PROPERTIES PUBLIC_HEADER ${analysis_SOURCE_DIR}/Lib/inc)

target_link_libraries(Interpolation ${ROOT_LIBRARIES})
#target_link_libraries(Selector ${ROOT_LIBRARIES} NPMugast NPMUST2 NPCATS)
target_link_libraries(Selector ${ROOT_LIBRARIES}  NPMugast NPMUST2 NPCATS Interpolation)
target_link_libraries(Analysis ${ROOT_LIBRARIES} Selector)

#install(FILES ${CMAKE_CURRENT_BINARY_DIR}/src/libInterpolation_rdict.pcm ${analysis_BUILD_DIR}/Lib/src/libInterpolation.rootmap DESTINATION ${analysis_BUILD_DIR}/Lib/bin/)
#install(FILES ${analysis_BUILD_DIR}/Lib/src/libSelector_rdict.pcm ${analysis_BUILD_DIR}/Lib/src/libSelector.rootmap DESTINATION ${analysis_BUILD_DIR}/Lib/bin/)
#install(FILES ${analysis_BUILD_DIR}/Lib/src/libAnalysis_rdict.pcm ${analysis_BUILD_DIR}/Lib/src/libAnalysis.rootmap DESTINATION ${analysis_BUILD_DIR}/Lib/bin/)

