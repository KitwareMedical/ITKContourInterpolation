itk_module_test()

set( MorphologicalContourInterpolationTests
  itkMorphologicalContourInterpolationTest.cxx
  )

CreateTestDriver( MorphologicalContourInterpolation "${MorphologicalContourInterpolation-Test_LIBRARIES}" "${MorphologicalContourInterpolationTests}" )

add_executable(runFromIDE manualTest.cxx ${MorphologicalContourInterpolationTests})
target_link_libraries(runFromIDE ${MorphologicalContourInterpolation-Test_LIBRARIES})

# used to produce an output which is then taken as a baseline for a regression test
function(JustInputTest ImageName Ext) #optional: Algo Axis Label
  foreach(f ${ARGN}) #create unique test name
    set(optArgs "${optArgs}_${f}")
  endforeach()
  set(outImage "${ITK_TEST_OUTPUT_DIR}/${ImageName}${optArgs}.${Ext}")
  # set(imageNameOpt "${ImageName}${optArgs}") # only used for regression
  
  itk_add_test( NAME itkMCI_JustInOut_${ImageName}${optArgs}
    COMMAND MorphologicalContourInterpolationTestDriver
    # --compare DATA{Baseline/${imageNameOpt}.${Ext}} ${outImage} # no regression
    itkMorphologicalContourInterpolationTest
      DATA{Input/${ImageName}.${Ext}} ${outImage} ${ARGN})
endfunction()

function(RegressionTest ImageName Ext) #optional: Algo Axis Label
  foreach(f ${ARGN}) #create unique test name
    set(optArgs "${optArgs}_${f}")
  endforeach()
  set(outImage "${ITK_TEST_OUTPUT_DIR}/${ImageName}${optArgs}.${Ext}")
  set(imageNameOpt "${ImageName}${optArgs}")
  
  itk_add_test( NAME itkMCI_${ImageName}${optArgs}
    COMMAND MorphologicalContourInterpolationTestDriver
    --compare DATA{Baseline/${imageNameOpt}.${Ext}} ${outImage}
    itkMorphologicalContourInterpolationTest
      DATA{Input/${ImageName}.${Ext}} ${outImage} ${ARGN})
endfunction()

# Handcrafted tests
RegressionTest(Empty nrrd ${MCI_Algorithm})
RegressionTest(NoSlices nrrd ${MCI_Algorithm}) #interpolation should not be invoked (output==input)
RegressionTest(Micro1 nrrd ${MCI_Algorithm})
RegressionTest(SimplestOneToOne nrrd ${MCI_Algorithm})
RegressionTest(OneToOne nrrd ${MCI_Algorithm})
RegressionTest(OneToThree nrrd ${MCI_Algorithm})
RegressionTest(ExtrapolationAppearing nrrd ${MCI_Algorithm})
RegressionTest(DoubleTwoLabelBranching nrrd ${MCI_Algorithm})
RegressionTest(TwoAxisDoubleTwoLabelBC nrrd ${MCI_Algorithm})
RegressionTest(AccidentalMiddleSliceSeg nrrd ${MCI_Algorithm})
RegressionTest(SevenLabels nrrd ${MCI_Algorithm})
RegressionTest(FullEnd nrrd ${MCI_Algorithm})
RegressionTest(1MN1 nrrd ${MCI_Algorithm}) # double sided many broken into many + 1
RegressionTest(ManyToMany16 nrrd ${MCI_Algorithm})
RegressionTest(ManyToMany nrrd ${MCI_Algorithm})
RegressionTest(GridSeg nrrd ${MCI_Algorithm})
RegressionTest(GridSeg2 nrrd ${MCI_Algorithm})

# Per axis/label tests
RegressionTest(GridSeg nrrd ${MCI_Algorithm} 0)
RegressionTest(GridSeg nrrd ${MCI_Algorithm} 1)
RegressionTest(GridSeg nrrd ${MCI_Algorithm} 2)
RegressionTest(GridSeg2 nrrd ${MCI_Algorithm} 0)
RegressionTest(GridSeg2 nrrd ${MCI_Algorithm} 1)
RegressionTest(GridSeg2 nrrd ${MCI_Algorithm} 2)
RegressionTest(ThreeAxisFourLabelConflict nrrd ${MCI_Algorithm})
RegressionTest(ThreeAxisFourLabelConflict nrrd ${MCI_Algorithm} 0)
RegressionTest(ThreeAxisFourLabelConflict nrrd ${MCI_Algorithm} 1)
RegressionTest(ThreeAxisFourLabelConflict nrrd ${MCI_Algorithm} 2)
RegressionTest(ThreeAxisFourLabelConflict nrrd ${MCI_Algorithm} 0 0) # equals to _0
RegressionTest(ThreeAxisFourLabelConflict nrrd ${MCI_Algorithm} 1 1)
RegressionTest(ThreeAxisFourLabelConflict nrrd ${MCI_Algorithm} 2 3) # equals to _2
RegressionTest(ThreeAxisFourLabelConflict nrrd ${MCI_Algorithm} 0 4) # equals to _0
RegressionTest(ThreeAxisFourLabelConflict nrrd ${MCI_Algorithm} 1 2)
RegressionTest(ThreeAxisFourLabelConflict nrrd ${MCI_Algorithm} 2 2) # does not modify image

# Paul's single axis tests
RegressionTest(105769moving nii.gz ${MCI_Algorithm}) #defaults: all axes all labels
RegressionTest(105769moving nii.gz ${MCI_Algorithm} 2) #axis 2
RegressionTest(105769moving nii.gz ${MCI_Algorithm} -1 0) #all axes all labels
RegressionTest(105769moving nii.gz ${MCI_Algorithm} -1 1) #all axes label 1
RegressionTest(105769moving nii.gz ${MCI_Algorithm} -1 2) #all axes label 2
RegressionTest(105769moving nii.gz ${MCI_Algorithm} 2 1) #axis 2 label 1
RegressionTest(105769moving nii.gz ${MCI_Algorithm} 2 2) #axis 2 label 2
RegressionTest(105769fixed nii.gz ${MCI_Algorithm}) #defaults: all axes all labels
RegressionTest(105769fixed nii.gz ${MCI_Algorithm} 2) #axis 2
RegressionTest(105769fixed nii.gz ${MCI_Algorithm} -1 0) #all axes all labels
RegressionTest(105769fixed nii.gz ${MCI_Algorithm} -1 1) #all axes label 1
RegressionTest(105769fixed nii.gz ${MCI_Algorithm} -1 2) #all axes label 2
RegressionTest(105769fixed nii.gz ${MCI_Algorithm} 2 1) #axis 2 label 1
RegressionTest(105769fixed nii.gz ${MCI_Algorithm} 2 2) #axis 2 label 2

# 64816L_amygdala_int
RegressionTest(64816L_amygdala_int nii.gz ${MCI_Algorithm})
RegressionTest(64816L_amygdala_int nii.gz ${MCI_Algorithm} 0)
RegressionTest(64816L_amygdala_int nii.gz ${MCI_Algorithm} 1)
RegressionTest(64816L_amygdala_int nii.gz ${MCI_Algorithm} 2)
RegressionTest(64816L_amygdala_int nii.gz ${MCI_Algorithm} 0 1)
RegressionTest(64816L_amygdala_int nii.gz ${MCI_Algorithm} 0 2)
RegressionTest(64816L_amygdala_int nii.gz ${MCI_Algorithm} 0 3)
RegressionTest(64816L_amygdala_int nii.gz ${MCI_Algorithm} 0 4)
RegressionTest(64816L_amygdala_int nii.gz ${MCI_Algorithm} 1 1)
RegressionTest(64816L_amygdala_int nii.gz ${MCI_Algorithm} 1 2)
RegressionTest(64816L_amygdala_int nii.gz ${MCI_Algorithm} 1 3)
RegressionTest(64816L_amygdala_int nii.gz ${MCI_Algorithm} 1 4)
RegressionTest(64816L_amygdala_int nii.gz ${MCI_Algorithm} 2 1)
RegressionTest(64816L_amygdala_int nii.gz ${MCI_Algorithm} 2 2)
RegressionTest(64816L_amygdala_int nii.gz ${MCI_Algorithm} 2 3)
RegressionTest(64816L_amygdala_int nii.gz ${MCI_Algorithm} 2 4)

# 64816L_amygdala_fin
RegressionTest(64816L_amygdala_fin nii.gz ${MCI_Algorithm})
RegressionTest(64816L_amygdala_fin nii.gz ${MCI_Algorithm} 0)
RegressionTest(64816L_amygdala_fin nii.gz ${MCI_Algorithm} 1)
RegressionTest(64816L_amygdala_fin nii.gz ${MCI_Algorithm} 2)
RegressionTest(64816L_amygdala_fin nii.gz ${MCI_Algorithm} 0 1)
RegressionTest(64816L_amygdala_fin nii.gz ${MCI_Algorithm} 0 2)
RegressionTest(64816L_amygdala_fin nii.gz ${MCI_Algorithm} 0 3)
RegressionTest(64816L_amygdala_fin nii.gz ${MCI_Algorithm} 0 4)
RegressionTest(64816L_amygdala_fin nii.gz ${MCI_Algorithm} 1 1)
RegressionTest(64816L_amygdala_fin nii.gz ${MCI_Algorithm} 1 2)
RegressionTest(64816L_amygdala_fin nii.gz ${MCI_Algorithm} 1 3)
RegressionTest(64816L_amygdala_fin nii.gz ${MCI_Algorithm} 1 4)
RegressionTest(64816L_amygdala_fin nii.gz ${MCI_Algorithm} 2 1)
RegressionTest(64816L_amygdala_fin nii.gz ${MCI_Algorithm} 2 2)
RegressionTest(64816L_amygdala_fin nii.gz ${MCI_Algorithm} 2 3)
RegressionTest(64816L_amygdala_fin nii.gz ${MCI_Algorithm} 2 4)

# 65239R_amygdala_int
RegressionTest(65239R_amygdala_int nii.gz ${MCI_Algorithm})
RegressionTest(65239R_amygdala_int nii.gz ${MCI_Algorithm} 0)
RegressionTest(65239R_amygdala_int nii.gz ${MCI_Algorithm} 1)
RegressionTest(65239R_amygdala_int nii.gz ${MCI_Algorithm} 2)
RegressionTest(65239R_amygdala_int nii.gz ${MCI_Algorithm} 0 1)
RegressionTest(65239R_amygdala_int nii.gz ${MCI_Algorithm} 0 2)
RegressionTest(65239R_amygdala_int nii.gz ${MCI_Algorithm} 0 3)
RegressionTest(65239R_amygdala_int nii.gz ${MCI_Algorithm} 0 4)
RegressionTest(65239R_amygdala_int nii.gz ${MCI_Algorithm} 1 1)
RegressionTest(65239R_amygdala_int nii.gz ${MCI_Algorithm} 1 2)
RegressionTest(65239R_amygdala_int nii.gz ${MCI_Algorithm} 1 3)
RegressionTest(65239R_amygdala_int nii.gz ${MCI_Algorithm} 1 4)
RegressionTest(65239R_amygdala_int nii.gz ${MCI_Algorithm} 2 1)
RegressionTest(65239R_amygdala_int nii.gz ${MCI_Algorithm} 2 2)
RegressionTest(65239R_amygdala_int nii.gz ${MCI_Algorithm} 2 3)
RegressionTest(65239R_amygdala_int nii.gz ${MCI_Algorithm} 2 4)

# 65239R_amygdala_fin
RegressionTest(65239R_amygdala_fin nii.gz ${MCI_Algorithm})
RegressionTest(65239R_amygdala_fin nii.gz ${MCI_Algorithm} 0)
RegressionTest(65239R_amygdala_fin nii.gz ${MCI_Algorithm} 1)
RegressionTest(65239R_amygdala_fin nii.gz ${MCI_Algorithm} 2)
RegressionTest(65239R_amygdala_fin nii.gz ${MCI_Algorithm} 0 1)
RegressionTest(65239R_amygdala_fin nii.gz ${MCI_Algorithm} 0 2)
RegressionTest(65239R_amygdala_fin nii.gz ${MCI_Algorithm} 0 3)
RegressionTest(65239R_amygdala_fin nii.gz ${MCI_Algorithm} 0 4)
RegressionTest(65239R_amygdala_fin nii.gz ${MCI_Algorithm} 1 1)
RegressionTest(65239R_amygdala_fin nii.gz ${MCI_Algorithm} 1 2)
RegressionTest(65239R_amygdala_fin nii.gz ${MCI_Algorithm} 1 3)
RegressionTest(65239R_amygdala_fin nii.gz ${MCI_Algorithm} 1 4)
RegressionTest(65239R_amygdala_fin nii.gz ${MCI_Algorithm} 2 1)
RegressionTest(65239R_amygdala_fin nii.gz ${MCI_Algorithm} 2 2)
RegressionTest(65239R_amygdala_fin nii.gz ${MCI_Algorithm} 2 3)
RegressionTest(65239R_amygdala_fin nii.gz ${MCI_Algorithm} 2 4)

# Long-running handcrafted tests
RegressionTest(BigZ nrrd ${MCI_Algorithm})
RegressionTest(BigZ nrrd ${MCI_Algorithm} 0)
RegressionTest(BigZ nrrd ${MCI_Algorithm} 1)
RegressionTest(BigZ nrrd ${MCI_Algorithm} 2)
RegressionTest(FaceRoI nrrd ${MCI_Algorithm})
RegressionTest(FaceRoI nrrd ${MCI_Algorithm} 0)
RegressionTest(FaceRoI nrrd ${MCI_Algorithm} 1)
RegressionTest(FaceRoI nrrd ${MCI_Algorithm} 2)

# Non-3D tests
RegressionTest(FourD nrrd ${MCI_Algorithm})
