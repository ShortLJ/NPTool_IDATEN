add_custom_command(OUTPUT NPLinearRansac3DDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh NPLinearRansac3D.h NPLinearRansac3DDict.cxx NPLinearRansac3D.rootmap libNPTrackReconstruction.so NPTrackReconstructionLinkDef.h DEPENDS NPLinearRansac3D.h)

add_custom_command(OUTPUT NPRansacDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh NPRansac.h NPRansacDict.cxx NPRansac.rootmap libNPTrackReconstruction.so NPTrackReconstructionLinkDef.h DEPENDS NPRansac.h)

add_custom_command(OUTPUT NPClusterDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh NPCluster.h NPClusterDict.cxx NPCluster.rootmap libNPTrackReconstruction.so NPTrackReconstructionLinkDef.h DEPENDS NPCluster.h)

add_custom_command(OUTPUT TrackingDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh Tracking.h TrackingDict.cxx Tracking.rootmap libNPTrackReconstruction.so NPTrackReconstructionLinkDef.h DEPENDS Tracking.h)

## Check if MINUIT2 is installed along with ROOT
find_library(libMinuit2_FOUND NAMES Minuit2 HINTS "${ROOTSYS}/lib")
if(libMinuit2_FOUND)
  message(STATUS "Minuit2 support enabled for TrackReconstruction.")
  add_definitions(-DHAVE_MINUIT2)
else()
  message(STATUS "Minuit2 support disabled for TrackReconstruction: *** NPDCReconstruction will not work at execution ***")
endif()      


add_library(NPTrackReconstruction SHARED NPTrackingUtility.cxx NPLinearRansac3D.cxx NPRansac.cxx NPCluster.cxx NPTrack.cxx Tracking.cxx NPLinearRansac3DDict.cxx NPRansacDict.cxx NPClusterDict.cxx TrackingDict.cxx NPDCReconstruction.cxx NPDCReconstructionMT.cxx)

if(Minuit2_FOUND)
target_link_libraries(NPTrackReconstruction ${ROOT_LIBRARIES} -lMinuit2 NPCore)
else()
target_link_libraries(NPTrackReconstruction ${ROOT_LIBRARIES} NPCore)
endif()

install(FILES NPLinearRansac3D.h NPRansac.h NPCluster.h NPTrack.h Tracking.h NPTrackingUtility.h NPDCReconstruction.h NPDCReconstructionMT.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})
