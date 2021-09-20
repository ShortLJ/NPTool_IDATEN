add_custom_command(OUTPUT TGAGGPhysicsDict.cxx COMMAND ../../scripts/build_dict.sh TGAGGPhysics.h TGAGGPhysicsDict.cxx TGAGGPhysics.rootmap libNPGAGG.dylib DEPENDS TGAGGPhysics.h)
add_custom_command(OUTPUT TGAGGDataDict.cxx COMMAND ../../scripts/build_dict.sh TGAGGData.h TGAGGDataDict.cxx TGAGGData.rootmap libNPGAGG.dylib DEPENDS TGAGGData.h)
add_library(NPGAGG SHARED TGAGGSpectra.cxx TGAGGData.cxx TGAGGPhysics.cxx TGAGGDataDict.cxx TGAGGPhysicsDict.cxx )
target_link_libraries(NPGAGG ${ROOT_LIBRARIES} NPCore) 
install(FILES TGAGGData.h TGAGGPhysics.h TGAGGSpectra.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

