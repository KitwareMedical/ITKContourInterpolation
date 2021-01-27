cmake_minimum_required(VERSION 3.10.2)

project(MorphologicalContourInterpolation_Examples CXX)

if(NOT ITK_SOURCE_DIR)
  find_package(ITK 5.2 REQUIRED COMPONENTS MorphologicalContourInterpolation ITKSmoothing ITKImageIO ITKTestKernel)
else()
  # when being built as part of ITK, ITKImageIO and ITKTransformIO
  # lists of modules are not yet ready, causing a configure error
  find_package(ITK REQUIRED COMPONENTS MorphologicalContourInterpolation ITKSmoothing ITKIOImageBase ITKIONRRD ITKTestKernel)
endif()
include(${ITK_USE_FILE})

add_executable(mciExample mciExample.cxx)
target_link_libraries(mciExample ${ITK_LIBRARIES})


# add some regression tests
enable_testing()

add_test(NAME mciExample_SevenLabels_3
  COMMAND mciExample ${CMAKE_CURRENT_SOURCE_DIR}/SevenLabels.nrrd
    ${CMAKE_CURRENT_BINARY_DIR}/SevenLabels_interpolated.nrrd
    3)

if(ITK_WRAP_PYTHON)
  add_test(NAME mciExample_ManyToMany
    COMMAND python ${CMAKE_CURRENT_LIST_DIR}/mciExample.py
      ${CMAKE_CURRENT_SOURCE_DIR}/ManyToMany.nrrd
      ${CMAKE_CURRENT_BINARY_DIR}/ManyToMany_interpolated.nrrd)
endif()