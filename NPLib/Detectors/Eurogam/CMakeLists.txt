add_custom_command(OUTPUT TEurogamDataDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TEurogamData.h TEurogamDataDict.cxx TEurogamData.rootmap libNPEurogam.dylib DEPENDS TEurogamData.h)
add_library(NPEurogam SHARED TEurogamData.cxx TEurogamDataDict.cxx )
target_link_libraries(NPEurogam ${ROOT_LIBRARIES} NPCore) 
install(FILES TEurogamData.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

