add_custom_command(OUTPUT TCharissaPhysicsDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TCharissaPhysics.h TCharissaPhysicsDict.cxx TCharissaPhysics.rootmap libNPCharissa.dylib DEPENDS TCharissaPhysics.h)
add_custom_command(OUTPUT TCharissaDataDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TCharissaData.h TCharissaDataDict.cxx TCharissaData.rootmap libNPCharissa.dylib DEPENDS TCharissaData.h)
add_library(NPCharissa SHARED TCharissaData.cxx TCharissaPhysics.cxx TCharissaDataDict.cxx TCharissaPhysicsDict.cxx TCharissaSpectra.cxx)
target_link_libraries(NPCharissa ${ROOT_LIBRARIES} NPCore) 
install(FILES TCharissaData.h TCharissaPhysics.h TCharissaSpectra.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

