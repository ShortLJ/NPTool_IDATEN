add_library(NPSFissionChamber SHARED  FissionChamber.cc)
target_link_libraries(NPSFissionChamber NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} ${NPLib_LIBRARIES} -lNPFissionChamber)
