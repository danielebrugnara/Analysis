ROOT_GENERATE_DICTIONARY(G__Plotter Plotter.h LINKDEF PlotterLinkDef.h)
add_library(Plotter STATIC Plotter.cxx G__Plotter.cxx)
target_link_libraries(Plotter ${ROOT_LIBRARIES} ReactionFragment ReactionReconstruction)

#linking to main
add_executable(kinematiclineplotter main.cxx)
target_link_libraries(kinematiclineplotter Plotter)
