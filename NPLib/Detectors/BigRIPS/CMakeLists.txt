add_custom_command(OUTPUT TBigRIPSPPACDataDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TBigRIPSPPACData.h TBigRIPSPPACDataDict.cxx TBigRIPSPPACData.rootmap libNPBigRIPS.dylib DEPENDS TBigRIPSPPACData.h)

add_custom_command(OUTPUT TBigRIPSPPACPhysicsDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TBigRIPSPPACPhysics.h TBigRIPSPPACPhysicsDict.cxx TBigRIPSPPACPhysics.rootmap libNPBigRIPS.dylib DEPENDS TBigRIPSPPACPhysics.h)

add_custom_command(OUTPUT TBigRIPSPlasticDataDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TBigRIPSPlasticData.h TBigRIPSPlasticDataDict.cxx TBigRIPSPlasticData.rootmap libNPBigRIPS.dylib DEPENDS TBigRIPSPlasticData.h)

add_custom_command(OUTPUT TBigRIPSPlasticPhysicsDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TBigRIPSPlasticPhysics.h TBigRIPSPlasticPhysicsDict.cxx TBigRIPSPlasticPhysics.rootmap libNPBigRIPS.dylib DEPENDS TBigRIPSPlasticPhysics.h)

add_custom_command(OUTPUT TBigRIPSICDataDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TBigRIPSICData.h TBigRIPSICDataDict.cxx TBigRIPSICData.rootmap libNPBigRIPS.dylib DEPENDS TBigRIPSICData.h)

add_custom_command(OUTPUT TBigRIPSICPhysicsDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TBigRIPSICPhysics.h TBigRIPSICPhysicsDict.cxx TBigRIPSICPhysics.rootmap libNPBigRIPS.dylib DEPENDS TBigRIPSICPhysics.h)

add_custom_command(OUTPUT TBigRIPSFocalDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TBigRIPSFocal.h TBigRIPSFocalDict.cxx TBigRIPSFocal.rootmap libNPBigRIPS.dylib DEPENDS TBigRIPSFocal.h)

add_custom_command(OUTPUT TBigRIPSRecoDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TBigRIPSReco.h TBigRIPSRecoDict.cxx TBigRIPSReco.rootmap libNPBigRIPS.dylib DEPENDS TBigRIPSReco.h)

add_library(NPBigRIPS SHARED TBigRIPSPPACData.cxx TBigRIPSPPACDataDict.cxx TBigRIPSPPACPhysics.cxx TBigRIPSPPACPhysicsDict.cxx TBigRIPSPlasticData.cxx TBigRIPSPlasticDataDict.cxx TBigRIPSPlasticPhysics.cxx TBigRIPSPlasticPhysicsDict.cxx TBigRIPSICData.cxx TBigRIPSICDataDict.cxx TBigRIPSICPhysics.cxx TBigRIPSICPhysicsDict.cxx TBigRIPSFocal.cxx TBigRIPSFocalDict.cxx TBigRIPSReco.cxx TBigRIPSRecoDict.cxx)

target_link_libraries(NPBigRIPS ${ROOT_LIBRARIES} NPCore NPTrackReconstruction) 
install(FILES TBigRIPSPPACData.h TBigRIPSPPACPhysics.h TBigRIPSPlasticData.h TBigRIPSPlasticPhysics.h TBigRIPSICData.h TBigRIPSICPhysics.h TBigRIPSFocal.h TBigRIPSReco.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

