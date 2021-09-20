add_library(NPSCatana SHARED  Catana.cc)
target_link_libraries(NPSCatana NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} ${NPLib_LIBRARIES} -lNPCatana)
