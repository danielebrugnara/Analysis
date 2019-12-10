add_subdirectory(src)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/src/libInterpolation_rdict.pcm ${CMAKE_CURRENT_BINARY_DIR}/src/libInterpolation.rootmap DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/bin/)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/src/libSelector_rdict.pcm ${CMAKE_CURRENT_BINARY_DIR}/src/libSelector.rootmap DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/bin/)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/src/libAnalysis_rdict.pcm ${CMAKE_CURRENT_BINARY_DIR}/src/libAnalysis.rootmap DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/bin/)