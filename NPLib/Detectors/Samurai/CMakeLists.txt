add_custom_command(OUTPUT TSamuraiHodoscopeDataDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TSamuraiHodoscopeData.h TSamuraiHodoscopeDataDict.cxx TSamuraiHodoscopeData.rootmap libNPSamurai.dylib DEPENDS TSamuraiHodoscopeData.h)


add_custom_command(OUTPUT SamuraiFieldMapDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh SamuraiFieldMap.h SamuraiFieldMapDict.cxx SamuraiFieldMap.rootmap libNPSamurai.dylib DEPENDS SamuraiFieldMap.h)

add_custom_command(OUTPUT TSamuraiHodoscopePhysicsDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TSamuraiHodoscopePhysics.h TSamuraiHodoscopePhysicsDict.cxx TSamuraiHodoscopePhysics.rootmap libNPSamurai.dylib DEPENDS TSamuraiHodoscopePhysics.h)

add_custom_command(OUTPUT TSamuraiBDCDataDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TSamuraiBDCData.h TSamuraiBDCDataDict.cxx TSamuraiBDCData.rootmap libNPSamurai.dylib DEPENDS TSamuraiBDCData.h)

add_custom_command(OUTPUT TSamuraiBDCPhysicsDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TSamuraiBDCPhysics.h TSamuraiBDCPhysicsDict.cxx TSamuraiBDCPhysics.rootmap libNPSamurai.dylib DEPENDS TSamuraiBDCPhysics.h)

add_custom_command(OUTPUT TSamuraiFDC2DataDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TSamuraiFDC2Data.h TSamuraiFDC2DataDict.cxx TSamuraiFDC2Data.rootmap libNPSamurai.dylib DEPENDS TSamuraiFDC2Data.h)

add_custom_command(OUTPUT TSamuraiFDC2PhysicsDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TSamuraiFDC2Physics.h TSamuraiFDC2PhysicsDict.cxx TSamuraiFDC2Physics.rootmap libNPSamurai.dylib DEPENDS TSamuraiFDC2Physics.h)

add_custom_command(OUTPUT TSamuraiFDC0DataDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TSamuraiFDC0Data.h TSamuraiFDC0DataDict.cxx TSamuraiFDC0Data.rootmap libNPSamurai.dylib DEPENDS TSamuraiFDC0Data.h)

add_custom_command(OUTPUT TSamuraiFDC0PhysicsDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TSamuraiFDC0Physics.h TSamuraiFDC0PhysicsDict.cxx TSamuraiFDC0Physics.rootmap libNPSamurai.dylib DEPENDS TSamuraiFDC0Physics.h)

add_library(NPSamurai SHARED SamuraiFieldMapDict.cxx SamuraiFieldMap.cxx TSamuraiHodoscopeData.cxx TSamuraiBDCData.cxx TSamuraiBDCPhysics.cxx TSamuraiBDCDataDict.cxx TSamuraiBDCPhysicsDict.cxx TSamuraiFDC2Data.cxx TSamuraiFDC2DataDict.cxx TSamuraiHodoscopePhysics.cxx TSamuraiHodoscopePhysicsDict.cxx TSamuraiFDC2Physics.cxx TSamuraiHodoscopeDataDict.cxx TSamuraiFDC2PhysicsDict.cxx TSamuraiFDC0Data.cxx TSamuraiFDC0DataDict.cxx TSamuraiFDC0Physics.cxx TSamuraiFDC0PhysicsDict.cxx)

target_link_libraries(NPSamurai ${ROOT_LIBRARIES} NPCore NPTrackReconstruction NPPhysics) 
install(FILES SamuraiFieldMap.h TSamuraiBDCData.h TSamuraiBDCPhysics.h TSamuraiHodoscopeData.h TSamuraiHodoscopePhysics.h TSamuraiFDC2Data.h TSamuraiFDC2Physics.h  TSamuraiFDC0Data.h TSamuraiFDC0Physics.h SamuraiDCIndex.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

