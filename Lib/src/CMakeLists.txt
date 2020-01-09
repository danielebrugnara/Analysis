include_directories(${analysis_SOURCE_DIR}/Lib/inc $ENV{NPTOOL}/NPLib/include)
#add_definitions(${ROOT_CXX_FLAGS})  

link_directories($ENV{NPTOOL}/NPLib/lib)


#Building root dictinaries
ROOT_GENERATE_DICTIONARY(G__Interpolation Interpolation.h LINKDEF ${analysis_SOURCE_DIR}/Lib/inc/InterpolationLinkDef.h)
add_library(Interpolation SHARED Interpolation.cxx G__Interpolation.cxx)

ROOT_GENERATE_DICTIONARY(G__Identification Identification.h LINKDEF ${analysis_SOURCE_DIR}/Lib/inc/IdentificationLinkDef.h)
ROOT_GENERATE_DICTIONARY(G__VamosIdentification VamosIdentification.h LINKDEF ${analysis_SOURCE_DIR}/Lib/inc/VamosIdentificationLinkDef.h)
add_library(VamosIdentification SHARED VamosIdentification.cxx Identification.cxx G__VamosIdentification.cxx G__Identification.cxx)

ROOT_GENERATE_DICTIONARY(G__Selector Selector.h LINKDEF ${analysis_SOURCE_DIR}/Lib/inc/SelectorLinkDef.h)
add_library(Selector SHARED Selector.cxx G__Selector.cxx)

ROOT_GENERATE_DICTIONARY(G__Analysis Analysis.h LINKDEF ${analysis_SOURCE_DIR}/Lib/inc/AnalysisLinkDef.h)
add_library(Analysis SHARED Analysis.cxx G__Analysis.cxx)


set_target_properties(Selector PROPERTIES PUBLIC_HEADER ${analysis_SOURCE_DIR}/Lib/inc)
set_target_properties(VamosIdentification PROPERTIES PUBLIC_HEADER ${analysis_SOURCE_DIR}/Lib/inc)
set_target_properties(Interpolation PROPERTIES PUBLIC_HEADER ${analysis_SOURCE_DIR}/Lib/inc)
set_target_properties(Analysis PROPERTIES PUBLIC_HEADER ${analysis_SOURCE_DIR}/Lib/inc)

#Linking to needed dynamic libraries
target_link_libraries(Interpolation ${ROOT_LIBRARIES})
target_link_libraries(VamosIdentification ${ROOT_LIBRARIES})
target_link_libraries(Selector ${ROOT_LIBRARIES}  NPMugast NPMUST2 NPCATS Interpolation VamosIdentification)
target_link_libraries(Analysis ${ROOT_LIBRARIES} Selector)
