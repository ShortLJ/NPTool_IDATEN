add_custom_command(OUTPUT TTACTICPhysicsDict.cxx COMMAND ../../scripts/build_dict.sh TTACTICPhysics.h TTACTICPhysicsDict.cxx TTACTICPhysics.rootmap libNPTACTIC.dylib DEPENDS TTACTICPhysics.h)
add_custom_command(OUTPUT TTACTICDataDict.cxx COMMAND ../../scripts/build_dict.sh TTACTICData.h TTACTICDataDict.cxx TTACTICData.rootmap libNPTACTIC.dylib DEPENDS TTACTICData.h)
add_library(NPTACTIC SHARED TTACTICSpectra.cxx TTACTICData.cxx TTACTICPhysics.cxx TTACTICDataDict.cxx TTACTICPhysicsDict.cxx )
target_link_libraries(NPTACTIC ${ROOT_LIBRARIES} NPCore) 
install(FILES TTACTICData.h TTACTICPhysics.h TTACTICSpectra.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

