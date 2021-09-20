add_library(NPSScone SHARED  Scone.cc)
target_link_libraries(NPSScone NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} ${NPLib_LIBRARIES} -lNPScone)
