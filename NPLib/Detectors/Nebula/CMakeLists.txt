add_custom_command(OUTPUT TNebulaPhysicsDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TNebulaPhysics.h TNebulaPhysicsDict.cxx TNebulaPhysics.rootmap libNPNebula.dylib DEPENDS TNebulaPhysics.h)
add_custom_command(OUTPUT TNebulaDataDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TNebulaData.h TNebulaDataDict.cxx TNebulaData.rootmap libNPNebula.dylib DEPENDS TNebulaData.h)
add_library(NPNebula SHARED TNebulaSpectra.cxx TNebulaData.cxx TNebulaPhysics.cxx TNebulaDataDict.cxx TNebulaPhysicsDict.cxx )
target_link_libraries(NPNebula ${ROOT_LIBRARIES} NPCore) 
install(FILES TNebulaData.h TNebulaPhysics.h TNebulaSpectra.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

