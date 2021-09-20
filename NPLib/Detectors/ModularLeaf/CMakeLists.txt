add_custom_command(OUTPUT TModularLeafPhysicsDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TModularLeafPhysics.h TModularLeafPhysicsDict.cxx TModularLeafPhysics.rootmap libNPModularLeaf.dylib DEPENDS TModularLeafPhysics.h)
add_library(NPModularLeaf SHARED TModularLeafPhysics.cxx TModularLeafPhysicsDict.cxx TModularLeafSpectra.cxx)
target_link_libraries(NPModularLeaf ${ROOT_LIBRARIES} NPCore) 
install(FILES TModularLeafPhysics.h TModularLeafSpectra.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

