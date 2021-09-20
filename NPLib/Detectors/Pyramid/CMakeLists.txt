add_custom_command(OUTPUT TPyramidPhysicsDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TPyramidPhysics.h TPyramidPhysicsDict.cxx TPyramidPhysics.rootmap libNPPyramid.dylib DEPENDS TPyramidPhysics.h)
add_custom_command(OUTPUT TPyramidDataDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TPyramidData.h TPyramidDataDict.cxx TPyramidData.rootmap libNPPyramid.dylib DEPENDS TPyramidData.h)
add_library(NPPyramid SHARED TPyramidSpectra.cxx TPyramidData.cxx TPyramidPhysics.cxx TPyramidDataDict.cxx TPyramidPhysicsDict.cxx )
target_link_libraries(NPPyramid ${ROOT_LIBRARIES} NPCore) 
install(FILES TPyramidData.h TPyramidPhysics.h TPyramidSpectra.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

