option(WITH_GARFIELD "Build with Garfield for electron drift in non-uniform field" OFF)

#IMPORTANT::IF CHANGING ABOVE OPTION DELETE CMakeCache IN nps FOR IT TO TAKE AFFECT. IF UPLOADING TO GIT SET THIS OPTION TO OFF, OTHERWISE A FULL COMPILATION OFF ALL DETECTORS IN NPTOOLS WILL REQUIRE GARFIELD!!! 

if(WITH_GARFIELD)
find_package(Garfield REQUIRED)
include_directories($ENV{GARFIELD_HOME}/Include)
link_directories($ENV{GARFIELD_HOME}/Library)
file(GLOB headers $ENV{GARFIELD_HOME}/Include/*.hh)
endif()

add_library(NPSTACTIC SHARED  TACTIC.cc TACTICScorer.cc)

if(WITH_GARFIELD)
target_link_libraries(NPSTACTIC NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} ${NPLib_LIBRARIES} -lNPTACTIC /Users/warren/garfieldpp/install/lib/libGarfield.0.3.0.dylib /Users/warren/garfieldpp/install/lib/libGarfield.dylib) 
else()
target_link_libraries(NPSTACTIC NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} ${NPLib_LIBRARIES} -lNPTACTIC)
endif()
