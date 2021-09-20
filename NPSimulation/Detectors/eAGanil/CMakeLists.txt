add_library(NPSeAGanil SHARED  eAGanil.cc)
target_link_libraries(NPSeAGanil NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} ${NPLib_LIBRARIES} -lNPeAGanil)
