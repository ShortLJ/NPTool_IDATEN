add_library(NPSGAGG SHARED  GAGG.cc)
target_link_libraries(NPSGAGG NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} ${NPLib_LIBRARIES} -lNPGAGG)
