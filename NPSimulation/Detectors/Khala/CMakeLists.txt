add_library(NPSKhala SHARED  Khala.cc)
target_link_libraries(NPSKhala NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} NPKhala)
