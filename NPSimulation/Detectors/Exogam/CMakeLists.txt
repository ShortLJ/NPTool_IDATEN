add_library(NPSExogam SHARED  Exogam.cc)
target_link_libraries(NPSExogam NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} ${NPLib_LIBRARIES} -lNPExogam)
