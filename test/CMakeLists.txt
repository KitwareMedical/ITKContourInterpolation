itk_module_test()

set( MorphologicalContourInterpolationTests
  itkMorphologicalContourInterpolationTest.cxx
  )

if(${MCI_TestWithRLEImage})
list(APPEND MorphologicalContourInterpolationTests
  itkMorphologicalContourInterpolationTestWithRLEImage.cxx
  )
endif(${MCI_TestWithRLEImage})

CreateTestDriver( MorphologicalContourInterpolation "${MorphologicalContourInterpolation-Test_LIBRARIES}" "${MorphologicalContourInterpolationTests}" )

add_executable(manualTest manualTest.cxx ${MorphologicalContourInterpolationTests})
target_link_libraries(manualTest ${MorphologicalContourInterpolation-Test_LIBRARIES})

add_executable(dscComparison dscComparison.cxx)
target_link_libraries(dscComparison ${MorphologicalContourInterpolation-Test_LIBRARIES})

# used to produce an output which is then taken as a baseline for a regression test
function(JustInputTest ImageName Ext) #optional: Axis Label
  foreach(f ${ARGN}) #create unique test name
    set(optArgs "${optArgs}_${f}")
  endforeach()
  foreach(algo B;C;T)
    set(outImage "${ITK_TEST_OUTPUT_DIR}/${ImageName}_${algo}${optArgs}.${Ext}")

    itk_add_test( NAME itkMCI_JustInOut_${ImageName}_${algo}${optArgs}
      COMMAND MorphologicalContourInterpolationTestDriver
      itkMorphologicalContourInterpolationTest
        DATA{Input/${ImageName}.${Ext}} ${outImage} ${algo} ${ARGN})
  endforeach()
endfunction()

function(RegressionTest RunsLong ImageName Ext) #optional: Axis Label
  foreach(f ${ARGN}) #create unique test name
    set(optArgs "${optArgs}_${f}")
  endforeach()

  # B = repeated dilations with ball structuring element
  # C = repeated dilations with cross structuring element
  # T = distance transform (not geodesic!)
  foreach(algo B;C;T)
    set(outImage "${ITK_TEST_OUTPUT_DIR}/${ImageName}_${algo}${optArgs}.${Ext}")
    set(imageNameOpt "${ImageName}_${algo}${optArgs}")

    if(${MCI_TestWithRLEImage})
      set(TestExe itkMorphologicalContourInterpolationTestWithRLEImage)
    else(${MCI_TestWithRLEImage})
      set(TestExe itkMorphologicalContourInterpolationTest)
    endif(${MCI_TestWithRLEImage})

    itk_add_test( NAME itkMCI_${ImageName}_${algo}${optArgs}
      COMMAND MorphologicalContourInterpolationTestDriver
      --compare DATA{Baseline/${imageNameOpt}.${Ext}} ${outImage}
      ${TestExe}
        DATA{Input/${ImageName}.${Ext}} ${outImage} ${algo} ${ARGN})
    if (RunsLong)
      set_property(TEST itkMCI_${ImageName}_${algo}${optArgs} APPEND PROPERTY LABELS RUNS_LONG)
    endif()
  endforeach()
endfunction()

# functional comparison using Dice Similarity Coefficient
function(DSCTest ImageName Ext)
  set(outFilenameBase "${ITK_TEST_OUTPUT_DIR}/${ImageName}")

  itk_add_test( NAME itkMCI_DSC_${ImageName}
    COMMAND dscComparison DATA{Input/${ImageName}.${Ext}} ${outFilenameBase} ${ARGN})
  set_property(TEST itkMCI_DSC_${ImageName} APPEND PROPERTY LABELS RUNS_LONG)
endfunction()

# Handcrafted tests
RegressionTest(0 Empty nrrd)
RegressionTest(0 NoSlices nrrd) #interpolation should not be invoked (output==input)
RegressionTest(0 Micro1 nrrd)
RegressionTest(0 SimplestOneToOne nrrd)
RegressionTest(0 OneToOne nrrd)
RegressionTest(0 OneToThree nrrd)
RegressionTest(0 ExtrapolationAppearing nrrd)
RegressionTest(0 DoubleTwoLabelBranching nrrd)
RegressionTest(0 TwoAxisDoubleTwoLabelBC nrrd)
RegressionTest(0 AccidentalMiddleSliceSeg nrrd)
RegressionTest(0 SevenLabels nrrd)
RegressionTest(0 FullEnd nrrd)
RegressionTest(0 1MN1 nrrd) # double sided many broken into many + 1
RegressionTest(0 ManyToMany16 nrrd)
RegressionTest(0 ManyToMany nrrd)
RegressionTest(0 GridSeg nrrd)
RegressionTest(0 GridSeg2 nrrd)

# Per axis/label tests
RegressionTest(0 GridSeg nrrd 0)
RegressionTest(0 GridSeg nrrd 1)
RegressionTest(0 GridSeg nrrd 2)
RegressionTest(0 GridSeg2 nrrd 0)
RegressionTest(0 GridSeg2 nrrd 1)
RegressionTest(0 GridSeg2 nrrd 2)
RegressionTest(0 ThreeAxisFourLabelConflict nrrd)
RegressionTest(0 ThreeAxisFourLabelConflict nrrd 0)
RegressionTest(0 ThreeAxisFourLabelConflict nrrd 1)
RegressionTest(0 ThreeAxisFourLabelConflict nrrd 2)
RegressionTest(0 ThreeAxisFourLabelConflict nrrd 0 0) # equals to _0
RegressionTest(0 ThreeAxisFourLabelConflict nrrd 1 1)
RegressionTest(0 ThreeAxisFourLabelConflict nrrd 2 3) # equals to _2
RegressionTest(0 ThreeAxisFourLabelConflict nrrd 0 4) # equals to _0
RegressionTest(0 ThreeAxisFourLabelConflict nrrd 1 2)
RegressionTest(0 ThreeAxisFourLabelConflict nrrd 2 2) # does not modify image

# Paul's single axis tests
RegressionTest(0 105769moving nii.gz) #defaults: all axes all labels
RegressionTest(0 105769moving nii.gz 2) #axis 2
RegressionTest(0 105769moving nii.gz -1 0) #all axes all labels
RegressionTest(0 105769moving nii.gz -1 1) #all axes label 1
RegressionTest(0 105769moving nii.gz -1 2) #all axes label 2
RegressionTest(0 105769moving nii.gz 2 1) #axis 2 label 1
RegressionTest(0 105769moving nii.gz 2 2) #axis 2 label 2
RegressionTest(0 105769fixed nii.gz) #defaults: all axes all labels
RegressionTest(0 105769fixed nii.gz 2) #axis 2
RegressionTest(0 105769fixed nii.gz -1 0) #all axes all labels
RegressionTest(0 105769fixed nii.gz -1 1) #all axes label 1
RegressionTest(0 105769fixed nii.gz -1 2) #all axes label 2
RegressionTest(0 105769fixed nii.gz 2 1) #axis 2 label 1
RegressionTest(0 105769fixed nii.gz 2 2) #axis 2 label 2

# 64816L_amygdala_int
RegressionTest(0 64816L_amygdala_int nii.gz)
RegressionTest(0 64816L_amygdala_int nii.gz 0)
RegressionTest(0 64816L_amygdala_int nii.gz 1)
RegressionTest(0 64816L_amygdala_int nii.gz 2)
RegressionTest(0 64816L_amygdala_int nii.gz 0 1)
RegressionTest(0 64816L_amygdala_int nii.gz 0 2)
RegressionTest(0 64816L_amygdala_int nii.gz 0 3)
RegressionTest(0 64816L_amygdala_int nii.gz 0 4)
RegressionTest(0 64816L_amygdala_int nii.gz 1 1)
RegressionTest(0 64816L_amygdala_int nii.gz 1 2)
RegressionTest(0 64816L_amygdala_int nii.gz 1 3)
RegressionTest(0 64816L_amygdala_int nii.gz 1 4)
RegressionTest(0 64816L_amygdala_int nii.gz 2 1)
RegressionTest(0 64816L_amygdala_int nii.gz 2 2)
RegressionTest(0 64816L_amygdala_int nii.gz 2 3)
RegressionTest(0 64816L_amygdala_int nii.gz 2 4)

# 64816L_amygdala_fin
RegressionTest(0 64816L_amygdala_fin nii.gz)
RegressionTest(0 64816L_amygdala_fin nii.gz 0)
RegressionTest(0 64816L_amygdala_fin nii.gz 1)
RegressionTest(0 64816L_amygdala_fin nii.gz 2)
RegressionTest(0 64816L_amygdala_fin nii.gz 0 1)
RegressionTest(0 64816L_amygdala_fin nii.gz 0 2)
RegressionTest(0 64816L_amygdala_fin nii.gz 0 3)
RegressionTest(0 64816L_amygdala_fin nii.gz 0 4)
RegressionTest(0 64816L_amygdala_fin nii.gz 1 1)
RegressionTest(0 64816L_amygdala_fin nii.gz 1 2)
RegressionTest(0 64816L_amygdala_fin nii.gz 1 3)
RegressionTest(0 64816L_amygdala_fin nii.gz 1 4)
RegressionTest(0 64816L_amygdala_fin nii.gz 2 1)
RegressionTest(0 64816L_amygdala_fin nii.gz 2 2)
RegressionTest(0 64816L_amygdala_fin nii.gz 2 3)
RegressionTest(0 64816L_amygdala_fin nii.gz 2 4)

# 65239R_amygdala_int
RegressionTest(0 65239R_amygdala_int nii.gz)
RegressionTest(0 65239R_amygdala_int nii.gz 0)
RegressionTest(0 65239R_amygdala_int nii.gz 1)
RegressionTest(0 65239R_amygdala_int nii.gz 2)
RegressionTest(0 65239R_amygdala_int nii.gz 0 1)
RegressionTest(0 65239R_amygdala_int nii.gz 0 2)
RegressionTest(0 65239R_amygdala_int nii.gz 0 3)
RegressionTest(0 65239R_amygdala_int nii.gz 0 4)
RegressionTest(0 65239R_amygdala_int nii.gz 1 1)
RegressionTest(0 65239R_amygdala_int nii.gz 1 2)
RegressionTest(0 65239R_amygdala_int nii.gz 1 3)
RegressionTest(0 65239R_amygdala_int nii.gz 1 4)
RegressionTest(0 65239R_amygdala_int nii.gz 2 1)
RegressionTest(0 65239R_amygdala_int nii.gz 2 2)
RegressionTest(0 65239R_amygdala_int nii.gz 2 3)
RegressionTest(0 65239R_amygdala_int nii.gz 2 4)

# 65239R_amygdala_fin
RegressionTest(0 65239R_amygdala_fin nii.gz)
RegressionTest(0 65239R_amygdala_fin nii.gz 0)
RegressionTest(0 65239R_amygdala_fin nii.gz 1)
RegressionTest(0 65239R_amygdala_fin nii.gz 2)
RegressionTest(0 65239R_amygdala_fin nii.gz 0 1)
RegressionTest(0 65239R_amygdala_fin nii.gz 0 2)
RegressionTest(0 65239R_amygdala_fin nii.gz 0 3)
RegressionTest(0 65239R_amygdala_fin nii.gz 0 4)
RegressionTest(0 65239R_amygdala_fin nii.gz 1 1)
RegressionTest(0 65239R_amygdala_fin nii.gz 1 2)
RegressionTest(0 65239R_amygdala_fin nii.gz 1 3)
RegressionTest(0 65239R_amygdala_fin nii.gz 1 4)
RegressionTest(0 65239R_amygdala_fin nii.gz 2 1)
RegressionTest(0 65239R_amygdala_fin nii.gz 2 2)
RegressionTest(0 65239R_amygdala_fin nii.gz 2 3)
RegressionTest(0 65239R_amygdala_fin nii.gz 2 4)

# Long-running handcrafted tests
RegressionTest(1 BigZ nrrd)
RegressionTest(1 BigZ nrrd 0)
RegressionTest(1 BigZ nrrd 1)
RegressionTest(1 BigZ nrrd 2)
RegressionTest(1 FaceRoI nrrd)
RegressionTest(1 FaceRoI nrrd 0)
RegressionTest(1 FaceRoI nrrd 1)
RegressionTest(1 FaceRoI nrrd 2)

# Non-3D tests
RegressionTest(0 FourD nrrd)

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
