
if(Geant4_VERSION_MAJOR GREATER 9)
  if(Geant4_VERSION_MINOR GREATER 4)

  message("Compiling with NeutronHPphysics")
  add_library(NPSProcess OBJECT FissionDecay.cc Decay.cc BeamReaction.cc FastDriftElectron.cc NPIonIonInelasticPhysic.cc PhysicsList.cc G4DriftElectron.cc G4IonizationWithDE.cc G4DriftElectronPhysics.cc G4DEAbsorption.cc G4DEAmplification.cc G4DETransport.cc menate_R.cc NeutronHPphysics.cc NeutronHPMessenger.cc)
    add_definitions(-DUSE_NEUTRONHP)
  else ()
  
  message("Compiling WITHOUT NeutronHPphysics, install geant4 10.5 or greater for this physics list")
  add_library(NPSProcess OBJECT FissionDecay.cc Decay.cc BeamReaction.cc FastDriftElectron.cc NPIonIonInelasticPhysic.cc PhysicsList.cc G4DriftElectron.cc G4IonizationWithDE.cc G4DriftElectronPhysics.cc G4DEAbsorption.cc G4DEAmplification.cc G4DETransport.cc menate_R.cc)
  endif()
 endif()


# target_link_libraries(NPSProcess ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} NPSCore ${NPLib_LIBRARIES} NPInitialConditions NPInteractionCoordinates NPReactionConditions NPFissionConditions)
