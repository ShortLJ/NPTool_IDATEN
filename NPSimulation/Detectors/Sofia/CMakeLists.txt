add_library(NPSSofia SHARED  Sofia.cc)
target_link_libraries(NPSSofia NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} ${NPLib_LIBRARIES} -lNPSofia)
