cmake_minimum_required(VERSION 3.16)
project(morphological-contour-interpolation)

set(CMAKE_CXX_STANDARD 20)

set(io_components)
if (NOT EMSCRIPTEN)
  set(io_components ITKIONIFTI)
endif()
find_package(ITK REQUIRED
  COMPONENTS
    WebAssemblyInterface
    MorphologicalContourInterpolation
    ${io_components}
    )
include(${ITK_USE_FILE})

add_executable(morphological-contour-interpolation morphological-contour-interpolation.cxx)
target_link_libraries(morphological-contour-interpolation PUBLIC ${ITK_LIBRARIES})

enable_testing()

add_test(NAME morphological-contour-interpolation COMMAND morphological-contour-interpolation
  ${CMAKE_CURRENT_SOURCE_DIR}/test/data/input/64816L_amygdala_int.nii.gz
  ${CMAKE_CURRENT_BINARY_DIR}/64816L_amygdala_int.nii.gz
  )