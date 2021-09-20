add_custom_command(OUTPUT TSweeperPhysicsDict.cxx COMMAND ../../scripts/build_dict.sh TSweeperPhysics.h TSweeperPhysicsDict.cxx TSweeperPhysics.rootmap libNPSweeper.dylib DEPENDS TSweeperPhysics.h)
add_custom_command(OUTPUT TSweeperDataDict.cxx COMMAND ../../scripts/build_dict.sh TSweeperData.h TSweeperDataDict.cxx TSweeperData.rootmap libNPSweeper.dylib DEPENDS TSweeperData.h)
add_library(NPSweeper SHARED TSweeperSpectra.cxx TSweeperData.cxx TSweeperPhysics.cxx TSweeperDataDict.cxx TSweeperPhysicsDict.cxx )
target_link_libraries(NPSweeper ${ROOT_LIBRARIES} NPCore) 
install(FILES TSweeperData.h TSweeperPhysics.h TSweeperSpectra.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

