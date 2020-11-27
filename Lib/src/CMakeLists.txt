include_directories(${analysis_SOURCE_DIR}/Lib/inc $ENV{NPTOOL}/NPLib/include)
#add_definitions(${ROOT_CXX_FLAGS})  

link_directories($ENV{NPTOOL}/NPLib/lib)

#Libraries without root
add_library(Units SHARED Units.cxx)
set_target_properties(Units PROPERTIES PUBLIC_HEADER ${analysis_SOURCE_DIR}/Lib/inc)

add_library(DiaGraph SHARED DiaGraph.cxx)
set_target_properties(DiaGraph PROPERTIES PUBLIC_HEADER ${analysis_SOURCE_DIR}/Lib/inc)

add_library(LevelScheme SHARED LevelScheme.cxx)
set_target_properties(LevelScheme PROPERTIES PUBLIC_HEADER ${analysis_SOURCE_DIR}/Lib/inc)
target_link_libraries(LevelScheme DiaGraph Units)


#Building root dictinaries
ROOT_GENERATE_DICTIONARY(G__Interpolation Interpolation.h LINKDEF ${analysis_SOURCE_DIR}/Lib/inc/InterpolationLinkDef.h)
add_library(Interpolation SHARED Interpolation.cxx G__Interpolation.cxx)
set_target_properties(Interpolation PROPERTIES PUBLIC_HEADER ${analysis_SOURCE_DIR}/Lib/inc)
target_link_libraries(Interpolation ${ROOT_LIBRARIES})

ROOT_GENERATE_DICTIONARY(G__Fitter Fitter.h LINKDEF ${analysis_SOURCE_DIR}/Lib/inc/FitterLinkDef.h)
add_library(Fitter SHARED Fitter.cxx G__Fitter.cxx)
set_target_properties(Fitter PROPERTIES PUBLIC_HEADER ${analysis_SOURCE_DIR}/Lib/inc)
target_link_libraries(Fitter ${ROOT_LIBRARIES})

ROOT_GENERATE_DICTIONARY(G__EnergyLoss EnergyLoss.h LINKDEF ${analysis_SOURCE_DIR}/Lib/inc/EnergyLossLinkDef.h)
add_library(EnergyLoss SHARED EnergyLoss.cxx G__EnergyLoss.cxx)
set_target_properties(EnergyLoss PROPERTIES PUBLIC_HEADER ${analysis_SOURCE_DIR}/Lib/inc)
target_link_libraries(EnergyLoss ${ROOT_LIBRARIES} Units)

ROOT_GENERATE_DICTIONARY(G__Minimizer Minimizer.h LINKDEF ${analysis_SOURCE_DIR}/Lib/inc/MinimizerLinkDef.h)
add_library(Minimizer SHARED Minimizer.cxx G__Minimizer.cxx)
set_target_properties(Minimizer PROPERTIES PUBLIC_HEADER ${analysis_SOURCE_DIR}/Lib/inc)

ROOT_GENERATE_DICTIONARY(G__ReactionFragment ReactionFragment.h LINKDEF ${analysis_SOURCE_DIR}/Lib/inc/ReactionFragmentLinkDef.h)
add_library(ReactionFragment SHARED ReactionFragment.cxx G__ReactionFragment.cxx)
set_target_properties(ReactionFragment PROPERTIES PUBLIC_HEADER ${analysis_SOURCE_DIR}/Lib/inc)
target_link_libraries(ReactionFragment ${ROOT_LIBRARIES})

ROOT_GENERATE_DICTIONARY(G__ReactionReconstruction ReactionReconstruction.h LINKDEF ${analysis_SOURCE_DIR}/Lib/inc/ReactionReconstructionLinkDef.h)
add_library(ReactionReconstruction SHARED ReactionReconstruction.cxx G__ReactionReconstruction.cxx)
set_target_properties(ReactionReconstruction PROPERTIES PUBLIC_HEADER ${analysis_SOURCE_DIR}/Lib/inc)
target_link_libraries(ReactionReconstruction ${ROOT_LIBRARIES} ReactionFragment)

ROOT_GENERATE_DICTIONARY(G__Calibration Calibration.h LINKDEF ${analysis_SOURCE_DIR}/Lib/inc/CalibrationLinkDef.h)
add_library(Calibration SHARED Calibration.cxx G__Calibration.cxx)
set_target_properties(Calibration PROPERTIES PUBLIC_HEADER ${analysis_SOURCE_DIR}/Lib/inc)
target_link_libraries(Calibration ${ROOT_LIBRARIES})

ROOT_GENERATE_DICTIONARY(G__VamosData VamosData.h LINKDEF ${analysis_SOURCE_DIR}/Lib/inc/VamosDataLinkDef.h)
add_library(VamosData SHARED VamosData.cxx G__VamosData.cxx)
set_target_properties(VamosData PROPERTIES PUBLIC_HEADER ${analysis_SOURCE_DIR}/Lib/inc)
target_link_libraries(VamosData ${ROOT_LIBRARIES})

ROOT_GENERATE_DICTIONARY(G__MugastData MugastData.h LINKDEF ${analysis_SOURCE_DIR}/Lib/inc/MugastDataLinkDef.h)
add_library(MugastData SHARED MugastData.cxx G__MugastData.cxx)
set_target_properties(MugastData PROPERTIES PUBLIC_HEADER ${analysis_SOURCE_DIR}/Lib/inc)
target_link_libraries(MugastData ${ROOT_LIBRARIES})

ROOT_GENERATE_DICTIONARY(G__Must2Data Must2Data.h LINKDEF ${analysis_SOURCE_DIR}/Lib/inc/Must2DataLinkDef.h)
add_library(Must2Data SHARED Must2Data.cxx G__Must2Data.cxx)
set_target_properties(Must2Data PROPERTIES PUBLIC_HEADER ${analysis_SOURCE_DIR}/Lib/inc)
target_link_libraries(Must2Data ${ROOT_LIBRARIES} MugastData)

ROOT_GENERATE_DICTIONARY(G__CatsData CatsData.h LINKDEF ${analysis_SOURCE_DIR}/Lib/inc/CatsDataLinkDef.h)
add_library(CatsData SHARED CatsData.cxx G__CatsData.cxx)
set_target_properties(CatsData PROPERTIES PUBLIC_HEADER ${analysis_SOURCE_DIR}/Lib/inc)
target_link_libraries(CatsData ${ROOT_LIBRARIES})

ROOT_GENERATE_DICTIONARY(G__AgataData AgataData.h LINKDEF ${analysis_SOURCE_DIR}/Lib/inc/AgataDataLinkDef.h)
add_library(AgataData SHARED AgataData.cxx G__AgataData.cxx)
set_target_properties(AgataData PROPERTIES PUBLIC_HEADER ${analysis_SOURCE_DIR}/Lib/inc)
target_link_libraries(AgataData ${ROOT_LIBRARIES})

ROOT_GENERATE_DICTIONARY(G__Identification Identification.h LINKDEF ${analysis_SOURCE_DIR}/Lib/inc/IdentificationLinkDef.h)

ROOT_GENERATE_DICTIONARY(G__VamosIdentification VamosIdentification.h LINKDEF ${analysis_SOURCE_DIR}/Lib/inc/VamosIdentificationLinkDef.h)
add_library(VamosIdentification SHARED VamosIdentification.cxx Identification.cxx G__VamosIdentification.cxx G__Identification.cxx)
set_target_properties(VamosIdentification PROPERTIES PUBLIC_HEADER ${analysis_SOURCE_DIR}/Lib/inc)
target_link_libraries(VamosIdentification ${ROOT_LIBRARIES} Interpolation VamosData)

ROOT_GENERATE_DICTIONARY(G__MugastIdentification MugastIdentification.h LINKDEF ${analysis_SOURCE_DIR}/Lib/inc/MugastIdentificationLinkDef.h)
add_library(MugastIdentification SHARED MugastIdentification.cxx Identification.cxx G__MugastIdentification.cxx G__Identification.cxx)
set_target_properties(MugastIdentification PROPERTIES PUBLIC_HEADER ${analysis_SOURCE_DIR}/Lib/inc)
target_link_libraries(MugastIdentification ${ROOT_LIBRARIES} Interpolation Minimizer NPMugast NPMUST2 NPCATS Calibration EnergyLoss MugastData Must2Data CatsData)

ROOT_GENERATE_DICTIONARY(G__AgataProcessing AgataProcessing.h LINKDEF ${analysis_SOURCE_DIR}/Lib/inc/AgataProcessingLinkDef.h)
add_library(AgataProcessing SHARED AgataProcessing.cxx G__AgataProcessing.cxx)
set_target_properties(AgataProcessing PROPERTIES PUBLIC_HEADER ${analysis_SOURCE_DIR}/Lib/inc)
target_link_libraries(AgataProcessing ${ROOT_LIBRARIES} AgataData)

ROOT_GENERATE_DICTIONARY(G__Selector Selector.h LINKDEF ${analysis_SOURCE_DIR}/Lib/inc/SelectorLinkDef.h)
add_library(Selector SHARED Selector.cxx G__Selector.cxx)
set_target_properties(Selector PROPERTIES PUBLIC_HEADER ${analysis_SOURCE_DIR}/Lib/inc)
target_link_libraries(Selector ${ROOT_LIBRARIES}  NPMugast NPMUST2 NPCATS VamosIdentification MugastIdentification AgataProcessing ReactionReconstruction)

ROOT_GENERATE_DICTIONARY(G__Analysis Analysis.h LINKDEF ${analysis_SOURCE_DIR}/Lib/inc/AnalysisLinkDef.h)
add_library(Analysis SHARED Analysis.cxx G__Analysis.cxx)
set_target_properties(Analysis PROPERTIES PUBLIC_HEADER ${analysis_SOURCE_DIR}/Lib/inc)
target_link_libraries(Analysis ${ROOT_LIBRARIES} Selector)

