add_library(NPSMugast SHARED  Mugast.cc)
target_link_libraries(NPSMugast NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} ${NPLib_LIBRARIES} -lNPMugast)
