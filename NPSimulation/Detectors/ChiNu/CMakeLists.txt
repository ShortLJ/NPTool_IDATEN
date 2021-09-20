add_library(NPSChiNu SHARED  ChiNu.cc)
target_link_libraries(NPSChiNu NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} ${NPLib_LIBRARIES} -lNPChiNu)
