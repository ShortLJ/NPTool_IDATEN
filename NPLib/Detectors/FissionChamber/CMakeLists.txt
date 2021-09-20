add_custom_command(OUTPUT TFissionChamberPhysicsDict.cxx COMMAND ../../scripts/build_dict.sh TFissionChamberPhysics.h TFissionChamberPhysicsDict.cxx TFissionChamberPhysics.rootmap libNPFissionChamber.dylib DEPENDS TFissionChamberPhysics.h)
add_custom_command(OUTPUT TFissionChamberDataDict.cxx COMMAND ../../scripts/build_dict.sh TFissionChamberData.h TFissionChamberDataDict.cxx TFissionChamberData.rootmap libNPFissionChamber.dylib DEPENDS TFissionChamberData.h)
add_library(NPFissionChamber SHARED TFissionChamberSpectra.cxx TFissionChamberData.cxx TFissionChamberPhysics.cxx TFissionChamberDataDict.cxx TFissionChamberPhysicsDict.cxx )
target_link_libraries(NPFissionChamber ${ROOT_LIBRARIES} NPCore) 
install(FILES TFissionChamberData.h TFissionChamberPhysics.h TFissionChamberSpectra.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

