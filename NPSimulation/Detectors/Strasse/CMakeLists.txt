add_library(NPSStrasse SHARED  Strasse.cc)
target_link_libraries(NPSStrasse NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} ${NPLib_LIBRARIES} -lNPStrasse)
