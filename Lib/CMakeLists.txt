add_subdirectory(src)

#Adding root dictionaries to the binary folder
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/src/libInterpolation_rdict.pcm ${CMAKE_CURRENT_BINARY_DIR}/src/libInterpolation.rootmap DESTINATION ${LIBRARY_OUTPUT_PATH})
##install(FILES ${CMAKE_CURRENT_BINARY_DIR}/src/libIdentification_rdict.pcm ${CMAKE_CURRENT_BINARY_DIR}/src/libIdentification.rootmap DESTINATION ${LIBRARY_OUTPUT_PATH})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/src/libVamosIdentification_rdict.pcm ${CMAKE_CURRENT_BINARY_DIR}/src/libVamosIdentification.rootmap DESTINATION ${LIBRARY_OUTPUT_PATH})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/src/libSelector_rdict.pcm ${CMAKE_CURRENT_BINARY_DIR}/src/libSelector.rootmap DESTINATION ${LIBRARY_OUTPUT_PATH})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/src/libAnalysis_rdict.pcm ${CMAKE_CURRENT_BINARY_DIR}/src/libAnalysis.rootmap DESTINATION ${LIBRARY_OUTPUT_PATH})
