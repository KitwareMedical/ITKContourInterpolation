cmake_minimum_required(VERSION 3.10.2)

project(MorphologicalContourInterpolation_Examples CXX)

find_package(ITK REQUIRED
  COMPONENTS
    MorphologicalContourInterpolation
    ITKSmoothing
    ITKIOImageBase
    ITKIONRRD
    )
include(${ITK_USE_FILE})

add_executable(mciExample mciExample.cxx)
target_link_libraries(mciExample ${ITK_LIBRARIES})


enable_testing()

add_test(NAME mciExample_SevenLabels_3
  COMMAND mciExample ${CMAKE_CURRENT_SOURCE_DIR}/SevenLabels.nrrd
    ${CMAKE_CURRENT_BINARY_DIR}/SevenLabels_interpolated.nrrd
    3)
