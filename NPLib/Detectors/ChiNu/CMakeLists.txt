add_custom_command(OUTPUT TChiNuPhysicsDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TChiNuPhysics.h TChiNuPhysicsDict.cxx TChiNuPhysics.rootmap libNPChiNu.dylib DEPENDS TChiNuPhysics.h)
add_custom_command(OUTPUT TChiNuDataDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TChiNuData.h TChiNuDataDict.cxx TChiNuData.rootmap libNPChiNu.dylib DEPENDS TChiNuData.h)
add_library(NPChiNu SHARED TChiNuSpectra.cxx TChiNuData.cxx TChiNuPhysics.cxx TChiNuDataDict.cxx TChiNuPhysicsDict.cxx )
target_link_libraries(NPChiNu ${ROOT_LIBRARIES} NPCore) 
install(FILES TChiNuData.h TChiNuPhysics.h TChiNuSpectra.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

