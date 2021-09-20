include(CheckLanguage)

add_custom_command(OUTPUT TMinosPhysicsDict.cxx COMMAND ../../scripts/build_dict.sh TMinosPhysics.h TMinosPhysicsDict.cxx TMinosPhysics.rootmap libNPMinos.dylib DEPENDS TMinosPhysics.h)
add_custom_command(OUTPUT TMinosDataDict.cxx COMMAND ../../scripts/build_dict.sh TMinosData.h TMinosDataDict.cxx TMinosData.rootmap libNPMinos.dylib DEPENDS TMinosData.h)

add_custom_command(OUTPUT oldTMinosDataDict.cxx COMMAND ../../scripts/build_dict.sh oldTMinosData.h oldTMinosDataDict.cxx oldTMinosData.rootmap libNPMinos.dylib DEPENDS oldTMinosData.h)

add_custom_command(OUTPUT oldTMinosPhysicsDict.cxx COMMAND ../../scripts/build_dict.sh oldTMinosPhysics.h oldTMinosPhysicsDict.cxx oldTMinosPhysics.rootmap libNPMinos.dylib DEPENDS oldTMinosPhysics.h)
add_custom_command(OUTPUT TMinosClustDict.cxx COMMAND ../../scripts/build_dict.sh TMinosClust.h TMinosClustDict.cxx TMinosClust.rootmap libNPMinos.dylib DEPENDS TMinosClust.h)
add_custom_command(OUTPUT TMinosResultDict.cxx COMMAND ../../scripts/build_dict.sh TMinosResult.h TMinosResultDict.cxx TMinosResult.rootmap libNPMinos.dylib DEPENDS TMinosResult.h)

## Check if MINUIT is installed along with ROOT
find_library(libMinuit_FOUND NAMES Minuit2 HINTS "$ENV{ROOTSYS}/lib")
if(libMinuit_FOUND)
  message(STATUS "Minuit support enabled for Minos.")
  add_definitions(-DHAVE_MINUIT)
else()
  message(STATUS "Minuit support disabled for Minos. Not found in $ENV{ROOTSYS}/lib")
endif()


if(libMinuit_FOUND)
## Link to Minuit library
  add_library(NPMinos SHARED TMinosClust.cxx TMinosData.cxx TMinosPhysics.cxx TMinosResult.cxx oldTMinosPhysicsDict.cxx oldTMinosDataDict.cxx TMinosDataDict.cxx TMinosPhysicsDict.cxx TMinosClustDict.cxx TMinosResultDict.cxx MinosUtility.cxx) 
  target_link_libraries(NPMinos ${ROOT_LIBRARIES} Minuit NPCore NPTrackReconstruction )
  install(FILES TMinosClust.h TMinosData.h TMinosPhysics.h TMinosResult.h MinosUtility.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})
else()
## No Minuit library, skip linking     
  add_library(NPMinos SHARED TMinosData.cxx TMinosDataDict.cxx) 
  target_link_libraries(NPMinos ${ROOT_LIBRARIES} NPCore )
  install(FILES TMinosData.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

endif()


