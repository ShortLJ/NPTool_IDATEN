add_custom_command(OUTPUT TStrassePhysicsDict.cxx COMMAND ../../scripts/build_dict.sh TStrassePhysics.h TStrassePhysicsDict.cxx TStrassePhysics.rootmap libNPStrasse.dylib DEPENDS TStrassePhysics.h)
add_custom_command(OUTPUT TStrasseDataDict.cxx COMMAND ../../scripts/build_dict.sh TStrasseData.h TStrasseDataDict.cxx TStrasseData.rootmap libNPStrasse.dylib DEPENDS TStrasseData.h)
add_library(NPStrasse SHARED TStrasseSpectra.cxx TStrasseData.cxx TStrassePhysics.cxx TStrasseDataDict.cxx TStrassePhysicsDict.cxx )
target_link_libraries(NPStrasse ${ROOT_LIBRARIES} NPCore) 
install(FILES TStrasseData.h TStrassePhysics.h TStrasseSpectra.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

