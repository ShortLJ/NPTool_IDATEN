add_custom_command(OUTPUT TPISTAPhysicsDict.cxx COMMAND ../../scripts/build_dict.sh TPISTAPhysics.h TPISTAPhysicsDict.cxx TPISTAPhysics.rootmap libNPPISTA.dylib DEPENDS TPISTAPhysics.h)

add_custom_command(OUTPUT TPISTADataDict.cxx COMMAND ../../scripts/build_dict.sh TPISTAData.h TPISTADataDict.cxx TPISTAData.rootmap libNPPISTA.dylib DEPENDS TPISTAData.h)

add_library(NPPISTA SHARED TPISTASpectra.cxx TPISTAData.cxx TPISTAPhysics.cxx TPISTADataDict.cxx TPISTAPhysicsDict.cxx )
target_link_libraries(NPPISTA ${ROOT_LIBRARIES} NPCore) 
install(FILES TPISTAData.h TPISTAPhysics.h TPISTASpectra.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

