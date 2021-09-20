add_custom_command(OUTPUT TKhalaPhysicsDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TKhalaPhysics.h TKhalaPhysicsDict.cxx TKhalaPhysics.rootmap libNPKhala.dylib DEPENDS TKhalaPhysics.h)
add_custom_command(OUTPUT TKhalaDataDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TKhalaData.h TKhalaDataDict.cxx TKhalaData.rootmap libNPKhala.dylib DEPENDS TKhalaData.h)
add_library(NPKhala SHARED TKhalaData.cxx TKhalaPhysics.cxx TKhalaDataDict.cxx TKhalaPhysicsDict.cxx )
target_link_libraries(NPKhala ${ROOT_LIBRARIES} NPCore) 
install(FILES TKhalaData.h TKhalaPhysics.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

