add_custom_command(OUTPUT TeAGanilPhysicsDict.cxx COMMAND ../../scripts/build_dict.sh TeAGanilPhysics.h TeAGanilPhysicsDict.cxx TeAGanilPhysics.rootmap libNPeAGanil.dylib DEPENDS TeAGanilPhysics.h)
add_custom_command(OUTPUT TeAGanilDataDict.cxx COMMAND ../../scripts/build_dict.sh TeAGanilData.h TeAGanilDataDict.cxx TeAGanilData.rootmap libNPeAGanil.dylib DEPENDS TeAGanilData.h)
add_library(NPeAGanil SHARED TeAGanilSpectra.cxx TeAGanilData.cxx TeAGanilPhysics.cxx TeAGanilDataDict.cxx TeAGanilPhysicsDict.cxx )
target_link_libraries(NPeAGanil ${ROOT_LIBRARIES} NPCore) 
install(FILES TeAGanilData.h TeAGanilPhysics.h TeAGanilSpectra.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

