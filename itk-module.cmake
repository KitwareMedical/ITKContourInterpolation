set(DOCUMENTATION "The modules provides filters to do interpolation
of manually segmented anatomical contours.
Enabling testing requires RLEImage module to be enabled.")

itk_module(MorphologicalContourInterpolation
  DEPENDS
    ITKBinaryMathematicalMorphology
    ITKDistanceMap
  TEST_DEPENDS
    ITKTestKernel
  EXCLUDE_FROM_DEFAULT
  DESCRIPTION
    "${DOCUMENTATION}"
  )
