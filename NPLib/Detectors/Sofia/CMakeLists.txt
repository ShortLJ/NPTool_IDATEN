add_custom_command(OUTPUT TSofSciDataDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TSofSciData.h TSofSciDataDict.cxx TSofSciData.rootmap libNPSofia.dylib DEPENDS TSofSciData.h)
add_custom_command(OUTPUT TSofSciPhysicsDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TSofSciPhysics.h TSofSciPhysicsDict.cxx TSofSciPhysics.rootmap libNPSofia.dylib DEPENDS TSofSciPhysics.h)

add_custom_command(OUTPUT TSofTofWDataDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TSofTofWData.h TSofTofWDataDict.cxx TSofTofWData.rootmap libNPSofia.dylib DEPENDS TSofTofWData.h)
add_custom_command(OUTPUT TSofTofWPhysicsDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TSofTofWPhysics.h TSofTofWPhysicsDict.cxx TSofTofWPhysics.rootmap libNPSofia.dylib DEPENDS TSofTofWPhysics.h)

add_custom_command(OUTPUT TSofMwpcDataDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TSofMwpcData.h TSofMwpcDataDict.cxx TSofMwpcData.rootmap libNPSofia.dylib DEPENDS TSofMwpcData.h)
add_custom_command(OUTPUT TSofMwpcPhysicsDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TSofMwpcPhysics.h TSofMwpcPhysicsDict.cxx TSofMwpcPhysics.rootmap libNPSofia.dylib DEPENDS TSofMwpcPhysics.h)

add_custom_command(OUTPUT TSofAtDataDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TSofAtData.h TSofAtDataDict.cxx TSofAtData.rootmap libNPSofia.dylib DEPENDS TSofAtData.h)
add_custom_command(OUTPUT TSofAtPhysicsDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TSofAtPhysics.h TSofAtPhysicsDict.cxx TSofAtPhysics.rootmap libNPSofia.dylib DEPENDS TSofAtPhysics.h)

add_custom_command(OUTPUT TSofTrimDataDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TSofTrimData.h TSofTrimDataDict.cxx TSofTrimData.rootmap libNPSofia.dylib DEPENDS TSofTrimData.h)
add_custom_command(OUTPUT TSofTrimPhysicsDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TSofTrimPhysics.h TSofTrimPhysicsDict.cxx TSofTrimPhysics.rootmap libNPSofia.dylib DEPENDS TSofTrimPhysics.h)

add_custom_command(OUTPUT TSofTwimDataDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TSofTwimData.h TSofTwimDataDict.cxx TSofTwimData.rootmap libNPSofia.dylib DEPENDS TSofTwimData.h)
add_custom_command(OUTPUT TSofTwimPhysicsDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TSofTwimPhysics.h TSofTwimPhysicsDict.cxx TSofTwimPhysics.rootmap libNPSofia.dylib DEPENDS TSofTwimPhysics.h)

add_custom_command(OUTPUT TSofBeamIDDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TSofBeamID.h TSofBeamIDDict.cxx TSofBeamID.rootmap libNPSofia.dylib DEPENDS TSofBeamID.h)

add_custom_command(OUTPUT TSofFissionFragmentDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TSofFissionFragment.h TSofFissionFragmentDict.cxx TSofFissionFragment.rootmap libNPSofia.dylib DEPENDS TSofFissionFragment.h)

add_library(NPSofia SHARED TSofSciData.cxx TSofSciDataDict.cxx TSofSciPhysics.cxx TSofSciPhysicsDict.cxx TSofMwpcData.cxx TSofMwpcDataDict.cxx TSofMwpcPhysics.cxx TSofMwpcPhysicsDict.cxx TSofAtData.cxx TSofAtDataDict.cxx TSofAtPhysics.cxx TSofAtPhysicsDict.cxx TSofTrimData.cxx TSofTrimDataDict.cxx TSofTrimPhysics.cxx TSofTrimPhysicsDict.cxx TSofTwimData.cxx TSofTwimDataDict.cxx TSofTwimPhysics.cxx TSofTwimPhysicsDict.cxx TSofTofWData.cxx TSofTofWDataDict.cxx TSofTofWPhysics.cxx TSofTofWPhysicsDict.cxx TSofBeamID.cxx TSofBeamIDDict.cxx TSofFissionFragment.cxx TSofFissionFragmentDict.cxx)

target_link_libraries(NPSofia ${ROOT_LIBRARIES} NPCore NPPhysics)

install(FILES TSofSciData.h TSofSciPhysics.h TSofMwpcData.h TSofMwpcPhysics.h TSofAtData.h TSofAtPhysics.h TSofTrimData.h TSofTrimPhysics.h TSofTwimData.h TSofTwimPhysics.h TSofTofWData.h TSofTofWPhysics.h TSofBeamID.h TSofFissionFragment.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

