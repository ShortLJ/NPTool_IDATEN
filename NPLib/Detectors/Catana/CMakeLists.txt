add_custom_command(OUTPUT TCatanaPhysicsDict.cxx COMMAND ../../scripts/build_dict.sh TCatanaPhysics.h TCatanaPhysicsDict.cxx TCatanaPhysics.rootmap libNPCatana.dylib DEPENDS TCatanaPhysics.h)
add_custom_command(OUTPUT TCatanaDataDict.cxx COMMAND ../../scripts/build_dict.sh TCatanaData.h TCatanaDataDict.cxx TCatanaData.rootmap libNPCatana.dylib DEPENDS TCatanaData.h)
add_library(NPCatana SHARED TCatanaSpectra.cxx TCatanaData.cxx TCatanaPhysics.cxx TCatanaDataDict.cxx TCatanaPhysicsDict.cxx )
target_link_libraries(NPCatana ${ROOT_LIBRARIES} NPCore NPTrackReconstruction) 
install(FILES TCatanaData.h TCatanaPhysics.h TCatanaSpectra.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

