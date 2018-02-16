cmake_minimum_required(VERSION 3.9.5)

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
find_path(INPUT_DATA_DIR
  NAMES ManyToMany16.nrrd
  HINTS ${ITK_DIR}/ExternalData/Modules/Remote/MorphologicalContourInterpolation/test/Input
        ${ITK_DIR}/ExternalData/Modules/External/MorphologicalContourInterpolation/test/Input
  )
set(TEST_OUTPUT_DIR ${MorphologicalContourInterpolation_Examples_BINARY_DIR}/Testing/Temporary)


if(INPUT_DATA_DIR)
  add_test(NAME mciExample_ManyToMany16_radius3
    COMMAND mciExample ${INPUT_DATA_DIR}/ManyToMany16.nrrd
      ${TEST_OUTPUT_DIR}/ManyToMany16_interpolated.nrrd
      3)
endif()