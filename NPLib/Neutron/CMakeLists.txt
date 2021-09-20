# add_custom_command(OUTPUT CrossTalkDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh CrossTalk.h CrossTalkDict.cxx CrossTalk.rootmap libNPNeutron.so NPNeutronLinkDef.h DEPENDS CrossTalk.h)

add_library(NPNeutron SHARED NPCrossTalk.cxx)

target_link_libraries(NPNeutron ${ROOT_LIBRARIES} NPCore)

install(FILES NPCrossTalk.h  DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})
