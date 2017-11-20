itk_module_test()

set( MorphologicalContourInterpolationTests
  itkMorphologicalContourInterpolationTest.cxx
  itkMorphologicalContourInterpolationTestWithRLEImage.cxx
  )

CreateTestDriver( MorphologicalContourInterpolation "${MorphologicalContourInterpolation-Test_LIBRARIES}" "${MorphologicalContourInterpolationTests}" )

add_executable(manualTest manualTest.cxx ${MorphologicalContourInterpolationTests})
target_link_libraries(manualTest ${MorphologicalContourInterpolation-Test_LIBRARIES})

add_executable(dscComparison dscComparison.cxx)
target_link_libraries(dscComparison ${MorphologicalContourInterpolation-Test_LIBRARIES})

# used to produce an output which is then taken as a baseline for a regression test
function(JustInputTest ImageName Ext) #optional: Axis Label
  foreach(f ${ARGN}) # for creating unique test name
    set(optArgs "${optArgs}_${f}")
  endforeach()
  foreach(algo B;C;T)
    set(imageAlgoArgs "${ImageName}_${algo}${optArgs}")

    itk_add_test( NAME itkMCI_JustInOut_${imageAlgoArgs}
      COMMAND MorphologicalContourInterpolationTestDriver
      itkMorphologicalContourInterpolationTest
      DATA{Input/${ImageName}.${Ext}} "${ITK_TEST_OUTPUT_DIR}/${imageAlgoArgs}.${Ext}" ${algo} ${ARGN})
  endforeach()
endfunction()

function(RegressionTest RunsLong ImageName Ext) #optional: Axis Label
  foreach(f ${ARGN}) # for creating unique test name
    set(optArgs "${optArgs}_${f}")
  endforeach()

  # B = repeated dilations with ball structuring element
  # C = repeated dilations with cross structuring element
  # T = distance transform (not geodesic!)
  foreach(algo B;C;T)
    set(imageAlgoArgs "${ImageName}_${algo}${optArgs}")
    set(outImage "${ITK_TEST_OUTPUT_DIR}/${imageAlgoArgs}.${Ext}")
    # have the RLE version output the results into a different file
    # otherwise parallel invocation of tests causes undeterministic failures
    # when both tests try to write to the same file at the same time
    set(outImageRLE "${ITK_TEST_OUTPUT_DIR}/${imageAlgoArgs}_RLE.${Ext}")

    itk_add_test( NAME itkMCI_${imageAlgoArgs}
      COMMAND MorphologicalContourInterpolationTestDriver
      --compare DATA{Baseline/${imageAlgoArgs}.${Ext}} ${outImage}
      itkMorphologicalContourInterpolationTest
      DATA{Input/${ImageName}.${Ext}} ${outImage} ${algo} ${ARGN})

    itk_add_test( NAME itkMCIRLE_${imageAlgoArgs}
      COMMAND MorphologicalContourInterpolationTestDriver
      --compare DATA{Baseline/${imageAlgoArgs}.${Ext}} ${outImageRLE}
      itkMorphologicalContourInterpolationTestWithRLEImage
      DATA{Input/${ImageName}.${Ext}} ${outImageRLE} ${algo} ${ARGN})

    set_property(TEST itkMCI_${imageAlgoArgs} APPEND PROPERTY LABELS MCI_PLAIN)
    set_property(TEST itkMCIRLE_${imageAlgoArgs} APPEND PROPERTY LABELS MCI_RLE)

    if (RunsLong)
      set_property(TEST itkMCI_${imageAlgoArgs} APPEND PROPERTY LABELS RUNS_LONG)
      set_property(TEST itkMCIRLE_${imageAlgoArgs} APPEND PROPERTY LABELS RUNS_LONG)
    endif()
  endforeach()
endfunction()

function(RegressionTestShort)
  RegressionTest(No ${ARGN})
endfunction()

function(RegressionTestLong)
  RegressionTest(Yes ${ARGN})
endfunction()

# functional comparison using Dice Similarity Coefficient
function(DSCTest ImageName Ext)
  set(outFilenameBase "${ITK_TEST_OUTPUT_DIR}/${ImageName}")

  itk_add_test( NAME itkMCI_DSC_${ImageName}
    COMMAND dscComparison DATA{Input/${ImageName}.${Ext}} ${outFilenameBase} ${ARGN})
  set_property(TEST itkMCI_DSC_${ImageName} APPEND PROPERTY LABELS RUNS_LONG)
endfunction()

# Handcrafted tests
RegressionTestShort(Empty nrrd)
RegressionTestShort(NoSlices nrrd) #interpolation should not be invoked (output==input)
RegressionTestShort(Micro1 nrrd)
RegressionTestShort(SimplestOneToOne nrrd)
RegressionTestShort(OneToOne nrrd)
RegressionTestShort(OneToThree nrrd)
RegressionTestShort(ExtrapolationAppearing nrrd)
RegressionTestShort(DoubleTwoLabelBranching nrrd)
RegressionTestShort(TwoAxisDoubleTwoLabelBC nrrd)
RegressionTestShort(AccidentalMiddleSliceSeg nrrd)
RegressionTestShort(SevenLabels nrrd)
RegressionTestShort(FullEnd nrrd)
RegressionTestShort(1MN1 nrrd) # double sided many broken into many + 1
RegressionTestShort(ManyToMany16 nrrd)
RegressionTestShort(ManyToMany nrrd)
RegressionTestShort(GridSeg nrrd)
RegressionTestShort(GridSeg2 nrrd)

# Per axis/label tests
RegressionTestShort(GridSeg nrrd 0)
RegressionTestShort(GridSeg nrrd 1)
RegressionTestShort(GridSeg nrrd 2)
RegressionTestShort(GridSeg2 nrrd 0)
RegressionTestShort(GridSeg2 nrrd 1)
RegressionTestShort(GridSeg2 nrrd 2)
RegressionTestShort(ThreeAxisFourLabelConflict nrrd)
RegressionTestShort(ThreeAxisFourLabelConflict nrrd 0)
RegressionTestShort(ThreeAxisFourLabelConflict nrrd 1)
RegressionTestShort(ThreeAxisFourLabelConflict nrrd 2)
RegressionTestShort(ThreeAxisFourLabelConflict nrrd 0 0) # equals to _0
RegressionTestShort(ThreeAxisFourLabelConflict nrrd 1 1)
RegressionTestShort(ThreeAxisFourLabelConflict nrrd 2 3) # equals to _2
RegressionTestShort(ThreeAxisFourLabelConflict nrrd 0 4) # equals to _0
RegressionTestShort(ThreeAxisFourLabelConflict nrrd 1 2)
RegressionTestShort(ThreeAxisFourLabelConflict nrrd 2 2) # does not modify image

# Paul's single axis tests
RegressionTestShort(105769moving nii.gz) #defaults: all axes all labels
RegressionTestShort(105769moving nii.gz 2) #axis 2
RegressionTestShort(105769moving nii.gz -1 0) #all axes all labels
RegressionTestShort(105769moving nii.gz -1 1) #all axes label 1
RegressionTestShort(105769moving nii.gz -1 2) #all axes label 2
RegressionTestShort(105769moving nii.gz 2 1) #axis 2 label 1
RegressionTestShort(105769moving nii.gz 2 2) #axis 2 label 2
RegressionTestShort(105769fixed nii.gz) #defaults: all axes all labels
RegressionTestShort(105769fixed nii.gz 2) #axis 2
RegressionTestShort(105769fixed nii.gz -1 0) #all axes all labels
RegressionTestShort(105769fixed nii.gz -1 1) #all axes label 1
RegressionTestShort(105769fixed nii.gz -1 2) #all axes label 2
RegressionTestShort(105769fixed nii.gz 2 1) #axis 2 label 1
RegressionTestShort(105769fixed nii.gz 2 2) #axis 2 label 2

# 64816L_amygdala_int
RegressionTestShort(64816L_amygdala_int nii.gz)
RegressionTestShort(64816L_amygdala_int nii.gz 0)
RegressionTestShort(64816L_amygdala_int nii.gz 1)
RegressionTestShort(64816L_amygdala_int nii.gz 2)
RegressionTestShort(64816L_amygdala_int nii.gz 0 1)
RegressionTestShort(64816L_amygdala_int nii.gz 0 2)
RegressionTestShort(64816L_amygdala_int nii.gz 0 3)
RegressionTestShort(64816L_amygdala_int nii.gz 0 4)
RegressionTestShort(64816L_amygdala_int nii.gz 1 1)
RegressionTestShort(64816L_amygdala_int nii.gz 1 2)
RegressionTestShort(64816L_amygdala_int nii.gz 1 3)
RegressionTestShort(64816L_amygdala_int nii.gz 1 4)
RegressionTestShort(64816L_amygdala_int nii.gz 2 1)
RegressionTestShort(64816L_amygdala_int nii.gz 2 2)
RegressionTestShort(64816L_amygdala_int nii.gz 2 3)
RegressionTestShort(64816L_amygdala_int nii.gz 2 4)

# 64816L_amygdala_fin
RegressionTestShort(64816L_amygdala_fin nii.gz)
RegressionTestShort(64816L_amygdala_fin nii.gz 0)
RegressionTestShort(64816L_amygdala_fin nii.gz 1)
RegressionTestShort(64816L_amygdala_fin nii.gz 2)
RegressionTestShort(64816L_amygdala_fin nii.gz 0 1)
RegressionTestShort(64816L_amygdala_fin nii.gz 0 2)
RegressionTestShort(64816L_amygdala_fin nii.gz 0 3)
RegressionTestShort(64816L_amygdala_fin nii.gz 0 4)
RegressionTestShort(64816L_amygdala_fin nii.gz 1 1)
RegressionTestShort(64816L_amygdala_fin nii.gz 1 2)
RegressionTestShort(64816L_amygdala_fin nii.gz 1 3)
RegressionTestShort(64816L_amygdala_fin nii.gz 1 4)
RegressionTestShort(64816L_amygdala_fin nii.gz 2 1)
RegressionTestShort(64816L_amygdala_fin nii.gz 2 2)
RegressionTestShort(64816L_amygdala_fin nii.gz 2 3)
RegressionTestShort(64816L_amygdala_fin nii.gz 2 4)

# 65239R_amygdala_int
RegressionTestShort(65239R_amygdala_int nii.gz)
RegressionTestShort(65239R_amygdala_int nii.gz 0)
RegressionTestShort(65239R_amygdala_int nii.gz 1)
RegressionTestShort(65239R_amygdala_int nii.gz 2)
RegressionTestShort(65239R_amygdala_int nii.gz 0 1)
RegressionTestShort(65239R_amygdala_int nii.gz 0 2)
RegressionTestShort(65239R_amygdala_int nii.gz 0 3)
RegressionTestShort(65239R_amygdala_int nii.gz 0 4)
RegressionTestShort(65239R_amygdala_int nii.gz 1 1)
RegressionTestShort(65239R_amygdala_int nii.gz 1 2)
RegressionTestShort(65239R_amygdala_int nii.gz 1 3)
RegressionTestShort(65239R_amygdala_int nii.gz 1 4)
RegressionTestShort(65239R_amygdala_int nii.gz 2 1)
RegressionTestShort(65239R_amygdala_int nii.gz 2 2)
RegressionTestShort(65239R_amygdala_int nii.gz 2 3)
RegressionTestShort(65239R_amygdala_int nii.gz 2 4)

# 65239R_amygdala_fin
RegressionTestShort(65239R_amygdala_fin nii.gz)
RegressionTestShort(65239R_amygdala_fin nii.gz 0)
RegressionTestShort(65239R_amygdala_fin nii.gz 1)
RegressionTestShort(65239R_amygdala_fin nii.gz 2)
RegressionTestShort(65239R_amygdala_fin nii.gz 0 1)
RegressionTestShort(65239R_amygdala_fin nii.gz 0 2)
RegressionTestShort(65239R_amygdala_fin nii.gz 0 3)
RegressionTestShort(65239R_amygdala_fin nii.gz 0 4)
RegressionTestShort(65239R_amygdala_fin nii.gz 1 1)
RegressionTestShort(65239R_amygdala_fin nii.gz 1 2)
RegressionTestShort(65239R_amygdala_fin nii.gz 1 3)
RegressionTestShort(65239R_amygdala_fin nii.gz 1 4)
RegressionTestShort(65239R_amygdala_fin nii.gz 2 1)
RegressionTestShort(65239R_amygdala_fin nii.gz 2 2)
RegressionTestShort(65239R_amygdala_fin nii.gz 2 3)
RegressionTestShort(65239R_amygdala_fin nii.gz 2 4)

# Long-running handcrafted tests
RegressionTestLong(BigZ nrrd)
RegressionTestLong(BigZ nrrd 0)
RegressionTestLong(BigZ nrrd 1)
RegressionTestLong(BigZ nrrd 2)
RegressionTestLong(FaceRoI nrrd)
RegressionTestLong(FaceRoI nrrd 0)
RegressionTestLong(FaceRoI nrrd 1)
RegressionTestLong(FaceRoI nrrd 2)

# Non-3D tests
RegressionTestShort(FourD nrrd)

# dscComparison (crops)
DSCTest(c2l  mha saveImages)
DSCTest(c2b  mha saveImages)
DSCTest(c10l mha saveImages)
DSCTest(c10b mha saveImages)

# dscComparison (full resolution)
DSCTest(case_2_labels  mha saveImages)
DSCTest(case_2_binary  mha saveImages)
DSCTest(case_10_labels mha saveImages)
DSCTest(case_10_binary mha saveImages)
