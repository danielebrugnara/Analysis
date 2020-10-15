#static library for selector

ROOT_GENERATE_DICTIONARY(G__Hit Hit.h LINKDEF HitLinkDef.h)
add_library(Hit SHARED  G__Hit.cxx)
target_link_libraries(Hit ${ROOT_LIBRARIES})

ROOT_GENERATE_DICTIONARY(G__Particle Particle.h LINKDEF ParticleLinkDef.h)
add_library(Particle SHARED  G__Particle.cxx)
target_link_libraries(Particle ${ROOT_LIBRARIES})

ROOT_GENERATE_DICTIONARY(G__Readsimu Readsimu.h LINKDEF ReadsimuLinkDef.h)
add_library(Readsimu SHARED G__Readsimu.cxx)
target_link_libraries(Readsimu ${ROOT_LIBRARIES} Hit Particle)

ROOT_GENERATE_DICTIONARY(G__Selector Selector.h LINKDEF SelectorLinkDef.h)
add_library(Selector STATIC Selector.cxx G__Selector.cxx)
target_link_libraries(Selector ${ROOT_LIBRARIES} Hit Particle DiaGraph)

ROOT_GENERATE_DICTIONARY(G__RunSelector RunSelector.h LINKDEF RunSelectorLinkDef.h)
add_library(RunSelector STATIC RunSelector.cxx G__RunSelector.cxx)
target_link_libraries(RunSelector ${ROOT_LIBRARIES} Fitter Selector Readsimu)

#linking to main
ROOT_GENERATE_DICTIONARY(G__main main.cxx LINKDEF mainLinkDef.h)
add_executable(simulationanalysis main.cxx)
target_link_libraries(simulationanalysis RunSelector ${ROOT_LIBRARIES})
