add_custom_command(OUTPUT NPParticleDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh NPParticle.h NPParticleDict.cxx NPParticle.rootmap libNPPhysics.so NPParticleLinkDef.h DEPENDS NPParticle.h NPParticleLinkDef.h) 

add_custom_command(OUTPUT NPReactionDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh NPReaction.h NPReactionDict.cxx NPReaction.rootmap libNPPhysics.so NPReactionLinkDef.h DEPENDS NPReaction.h NPReactionLinkDef.h) 

add_custom_command(OUTPUT NPQFSDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh NPQFS.h NPQFSDict.cxx NPQFS.rootmap libNPPhysics.so NPQFSLinkDef.h DEPENDS NPQFS.h NPQFSLinkDef.h) 

add_custom_command(OUTPUT NPDecayDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh NPDecay.h NPDecayDict.cxx NPDecay.rootmap libNPPhysics.so NPDecayLinkDef.h DEPENDS NPDecay.h NPDecayLinkDef.h) 

add_custom_command(OUTPUT NPFissionDecayDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh NPFissionDecay.h NPFissionDecayDict.cxx NPFissionDecay.rootmap libNPPhysics.so NPFissionDecayLinkDef.h DEPENDS NPFissionDecay.h NPFissionDecayLinkDef.h) 

add_custom_command(OUTPUT GEFDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh GEF.h GEFDict.cxx GEF.rootmap libNPPhysics.so GEFLinkDef.h DEPENDS GEF.h GEFLinkDef.h) 

add_custom_command(OUTPUT NPEnergyLossDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh NPEnergyLoss.h NPEnergyLossDict.cxx NPEnergyLoss.rootmap libNPPhysics.so NPEnergyLossLinkDef.h DEPENDS NPEnergyLoss.h NPEnergyLossLinkDef.h) 

add_custom_command(OUTPUT TInitialConditionsDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TInitialConditions.h TInitialConditionsDict.cxx TInitialConditions.rootmap libNPInitialConditions.so DEPENDS TInitialConditions.h)

add_custom_command(OUTPUT TTrackInfoDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TTrackInfo.h TTrackInfoDict.cxx TTrackInfo.rootmap libNPTrackInfo.so DEPENDS TTrackInfo.h)

add_custom_command(OUTPUT TInteractionCoordinatesDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TInteractionCoordinates.h TInteractionCoordinatesDict.cxx TInteractionCoordinates.rootmap libNPInteractionCoordinates.so DEPENDS TInteractionCoordinates.h)

add_custom_command(OUTPUT TReactionConditionsDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TReactionConditions.h TReactionConditionsDict.cxx TReactionConditions.rootmap libNPReactionConditions.so DEPENDS TReactionConditions.h)

add_custom_command(OUTPUT TFissionConditionsDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TFissionConditions.h TFissionConditionsDict.cxx TFissionConditions.rootmap libNPFissionConditions.so DEPENDS TFissionConditions.h)

add_library(NPPhysics SHARED GEF.cxx NPFissionDecay.cxx NPDecay.cxx NPBeam.cxx NPEnergyLoss.cxx NPFunction.cxx NPParticle.cxx NPReaction.cxx NPQFS.cxx NPParticleDict.cxx NPReactionDict.cxx NPQFSDict.cxx NPEnergyLossDict.cxx )
target_link_libraries(NPPhysics ${ROOT_LIBRARIES} Physics NPCore) 

add_library(NPInitialConditions  SHARED  TInitialConditions.cxx TInitialConditionsDict.cxx )
target_link_libraries(NPInitialConditions  ${ROOT_LIBRARIES} ) 

add_library(NPTrackInfo  SHARED  TTrackInfo.cxx TTrackInfoDict.cxx )
target_link_libraries(NPTrackInfo  ${ROOT_LIBRARIES} ) 

add_library(NPInteractionCoordinates SHARED TInteractionCoordinates.cxx TInteractionCoordinatesDict.cxx)
target_link_libraries(NPInteractionCoordinates ${ROOT_LIBRARIES} ) 

add_library(NPReactionConditions SHARED TReactionConditions.cxx TReactionConditionsDict.cxx)
target_link_libraries(NPReactionConditions ${ROOT_LIBRARIES} ) 

add_library(NPFissionConditions SHARED TFissionConditions.cxx TFissionConditionsDict.cxx)
target_link_libraries(NPFissionConditions ${ROOT_LIBRARIES} ) 


install(FILES GEF.h NPFissionDecay.h NPDecay.h NPBeam.h NPEnergyLoss.h NPFunction.h NPParticle.h NPNucleus.h NPReaction.h  NPQFS.h TInitialConditions.h TTrackInfo.h TInteractionCoordinates.h TReactionConditions.h TFissionConditions.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})
