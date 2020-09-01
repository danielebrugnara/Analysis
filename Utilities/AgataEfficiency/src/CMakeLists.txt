#static library for selector
ROOT_GENERATE_DICTIONARY(G__Globals Globals.h LINKDEF GlobalsLinkDef.h)
add_library(Globals STATIC G__Globals.cxx)
target_link_libraries(Globals ${ROOT_LIBRARIES})


ROOT_GENERATE_DICTIONARY(G__Plotter Plotter.h LINKDEF PlotterLinkDef.h)
add_library(Plotter STATIC Plotter.cxx G__Plotter.cxx)
target_link_libraries(Plotter ${ROOT_LIBRARIES})

ROOT_GENERATE_DICTIONARY(G__Fitter Fitter.h LINKDEF FitterLinkDef.h)
add_library(Fitter STATIC Fitter.cxx G__Fitter.cxx)
target_link_libraries(Fitter ${ROOT_LIBRARIES})

ROOT_GENERATE_DICTIONARY(G__Selector Selector.h LINKDEF SelectorLinkDef.h)
add_library(Selector STATIC Selector.cxx G__Selector.cxx)
target_link_libraries(Selector ${ROOT_LIBRARIES})

ROOT_GENERATE_DICTIONARY(G__RunSelector RunSelector.h LINKDEF RunSelectorLinkDef.h)
add_library(RunSelector STATIC RunSelector.cxx G__RunSelector.cxx)
target_link_libraries(RunSelector ${ROOT_LIBRARIES} Selector)


ROOT_GENERATE_DICTIONARY(G__SpectrumAnalyzer SpectrumAnalyzer.h LINKDEF SpectrumAnalyzerLinkDef.h)
add_library(SpectrumAnalyzer STATIC SpectrumAnalyzer.cxx G__SpectrumAnalyzer.cxx)
target_link_libraries(SpectrumAnalyzer ${ROOT_LIBRARIES} DiaGraph LevelScheme Plotter Fitter)

#linking to main
add_executable(effanalysis main.cxx)
target_link_libraries(effanalysis RunSelector SpectrumAnalyzer)
