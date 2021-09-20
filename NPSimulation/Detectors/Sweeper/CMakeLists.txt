add_library(NPSSweeper SHARED  Sweeper.cc SweeperScorers.cc MagField.cc)
target_link_libraries(NPSSweeper NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} ${NPLib_LIBRARIES} -lNPSweeper)
