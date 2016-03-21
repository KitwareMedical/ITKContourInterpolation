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
RegressionTest(Empty nrrd RD)
RegressionTest(NoSlices nrrd RD) #interpolation should not be invoked (output==input)
RegressionTest(Micro1 nrrd RD)
RegressionTest(SimplestOneToOne nrrd RD)
RegressionTest(OneToOne nrrd RD)
RegressionTest(OneToThree nrrd RD)
RegressionTest(ExtrapolationAppearing nrrd RD)
RegressionTest(DoubleTwoLabelBranching nrrd RD)
RegressionTest(TwoAxisDoubleTwoLabelBC nrrd RD)
RegressionTest(AccidentalMiddleSliceSeg nrrd RD)
RegressionTest(SevenLabels nrrd RD)
RegressionTest(FullEnd nrrd RD)
RegressionTest(1MN1 nrrd RD) # double sided many broken into many + 1
RegressionTest(ManyToMany16 nrrd RD)
RegressionTest(ManyToMany nrrd RD)
RegressionTest(GridSeg nrrd RD)
RegressionTest(GridSeg2 nrrd RD)

# Per axis/label tests
RegressionTest(GridSeg nrrd RD 0)
RegressionTest(GridSeg nrrd RD 1)
RegressionTest(GridSeg nrrd RD 2)
RegressionTest(GridSeg2 nrrd RD 0)
RegressionTest(GridSeg2 nrrd RD 1)
RegressionTest(GridSeg2 nrrd RD 2)
RegressionTest(ThreeAxisFourLabelConflict nrrd RD)
RegressionTest(ThreeAxisFourLabelConflict nrrd RD 0)
RegressionTest(ThreeAxisFourLabelConflict nrrd RD 1)
RegressionTest(ThreeAxisFourLabelConflict nrrd RD 2)
RegressionTest(ThreeAxisFourLabelConflict nrrd RD 0 0) # equals to _0
RegressionTest(ThreeAxisFourLabelConflict nrrd RD 1 1)
RegressionTest(ThreeAxisFourLabelConflict nrrd RD 2 3) # equals to _2
RegressionTest(ThreeAxisFourLabelConflict nrrd RD 0 4) # equals to _0
RegressionTest(ThreeAxisFourLabelConflict nrrd RD 1 2)
RegressionTest(ThreeAxisFourLabelConflict nrrd RD 2 2) # does not modify image

# Paul's single axis tests
RegressionTest(105769moving nii.gz RD) #defaults: all axes all labels
RegressionTest(105769moving nii.gz RD 2) #axis 2
RegressionTest(105769moving nii.gz RD -1 0) #all axes all labels
RegressionTest(105769moving nii.gz RD -1 1) #all axes label 1
RegressionTest(105769moving nii.gz RD -1 2) #all axes label 2
RegressionTest(105769moving nii.gz RD 2 1) #axis 2 label 1
RegressionTest(105769moving nii.gz RD 2 2) #axis 2 label 2
RegressionTest(105769fixed nii.gz RD) #defaults: all axes all labels
RegressionTest(105769fixed nii.gz RD 2) #axis 2
RegressionTest(105769fixed nii.gz RD -1 0) #all axes all labels
RegressionTest(105769fixed nii.gz RD -1 1) #all axes label 1
RegressionTest(105769fixed nii.gz RD -1 2) #all axes label 2
RegressionTest(105769fixed nii.gz RD 2 1) #axis 2 label 1
RegressionTest(105769fixed nii.gz RD 2 2) #axis 2 label 2

# 64816L_amygdala_int
RegressionTest(64816L_amygdala_int nii.gz RD)
RegressionTest(64816L_amygdala_int nii.gz RD 0)
RegressionTest(64816L_amygdala_int nii.gz RD 1)
RegressionTest(64816L_amygdala_int nii.gz RD 2)
RegressionTest(64816L_amygdala_int nii.gz RD 0 1)
RegressionTest(64816L_amygdala_int nii.gz RD 0 2)
RegressionTest(64816L_amygdala_int nii.gz RD 0 3)
RegressionTest(64816L_amygdala_int nii.gz RD 0 4)
RegressionTest(64816L_amygdala_int nii.gz RD 1 1)
RegressionTest(64816L_amygdala_int nii.gz RD 1 2)
RegressionTest(64816L_amygdala_int nii.gz RD 1 3)
RegressionTest(64816L_amygdala_int nii.gz RD 1 4)
RegressionTest(64816L_amygdala_int nii.gz RD 2 1)
RegressionTest(64816L_amygdala_int nii.gz RD 2 2)
RegressionTest(64816L_amygdala_int nii.gz RD 2 3)
RegressionTest(64816L_amygdala_int nii.gz RD 2 4)

# 64816L_amygdala_fin
RegressionTest(64816L_amygdala_fin nii.gz RD)
RegressionTest(64816L_amygdala_fin nii.gz RD 0)
RegressionTest(64816L_amygdala_fin nii.gz RD 1)
RegressionTest(64816L_amygdala_fin nii.gz RD 2)
RegressionTest(64816L_amygdala_fin nii.gz RD 0 1)
RegressionTest(64816L_amygdala_fin nii.gz RD 0 2)
RegressionTest(64816L_amygdala_fin nii.gz RD 0 3)
RegressionTest(64816L_amygdala_fin nii.gz RD 0 4)
RegressionTest(64816L_amygdala_fin nii.gz RD 1 1)
RegressionTest(64816L_amygdala_fin nii.gz RD 1 2)
RegressionTest(64816L_amygdala_fin nii.gz RD 1 3)
RegressionTest(64816L_amygdala_fin nii.gz RD 1 4)
RegressionTest(64816L_amygdala_fin nii.gz RD 2 1)
RegressionTest(64816L_amygdala_fin nii.gz RD 2 2)
RegressionTest(64816L_amygdala_fin nii.gz RD 2 3)
RegressionTest(64816L_amygdala_fin nii.gz RD 2 4)

# 65239R_amygdala_int
RegressionTest(65239R_amygdala_int nii.gz RD)
RegressionTest(65239R_amygdala_int nii.gz RD 0)
RegressionTest(65239R_amygdala_int nii.gz RD 1)
RegressionTest(65239R_amygdala_int nii.gz RD 2)
RegressionTest(65239R_amygdala_int nii.gz RD 0 1)
RegressionTest(65239R_amygdala_int nii.gz RD 0 2)
RegressionTest(65239R_amygdala_int nii.gz RD 0 3)
RegressionTest(65239R_amygdala_int nii.gz RD 0 4)
RegressionTest(65239R_amygdala_int nii.gz RD 1 1)
RegressionTest(65239R_amygdala_int nii.gz RD 1 2)
RegressionTest(65239R_amygdala_int nii.gz RD 1 3)
RegressionTest(65239R_amygdala_int nii.gz RD 1 4)
RegressionTest(65239R_amygdala_int nii.gz RD 2 1)
RegressionTest(65239R_amygdala_int nii.gz RD 2 2)
RegressionTest(65239R_amygdala_int nii.gz RD 2 3)
RegressionTest(65239R_amygdala_int nii.gz RD 2 4)

# 65239R_amygdala_fin
RegressionTest(65239R_amygdala_fin nii.gz RD)
RegressionTest(65239R_amygdala_fin nii.gz RD 0)
RegressionTest(65239R_amygdala_fin nii.gz RD 1)
RegressionTest(65239R_amygdala_fin nii.gz RD 2)
RegressionTest(65239R_amygdala_fin nii.gz RD 0 1)
RegressionTest(65239R_amygdala_fin nii.gz RD 0 2)
RegressionTest(65239R_amygdala_fin nii.gz RD 0 3)
RegressionTest(65239R_amygdala_fin nii.gz RD 0 4)
RegressionTest(65239R_amygdala_fin nii.gz RD 1 1)
RegressionTest(65239R_amygdala_fin nii.gz RD 1 2)
RegressionTest(65239R_amygdala_fin nii.gz RD 1 3)
RegressionTest(65239R_amygdala_fin nii.gz RD 1 4)
RegressionTest(65239R_amygdala_fin nii.gz RD 2 1)
RegressionTest(65239R_amygdala_fin nii.gz RD 2 2)
RegressionTest(65239R_amygdala_fin nii.gz RD 2 3)
RegressionTest(65239R_amygdala_fin nii.gz RD 2 4)

# Long-running handcrafted tests
RegressionTest(BigZ nrrd RD)
RegressionTest(BigZ nrrd RD 0)
RegressionTest(BigZ nrrd RD 1)
RegressionTest(BigZ nrrd RD 2)
RegressionTest(FaceRoI nrrd RD)
RegressionTest(FaceRoI nrrd RD 0)
RegressionTest(FaceRoI nrrd RD 1)
RegressionTest(FaceRoI nrrd RD 2)

# Non-3D tests
RegressionTest(TwoD png RD)
RegressionTest(TwoD png RD -1 255) # same as regular
RegressionTest(TwoD png RD 0) # same as regular
RegressionTest(TwoD png RD -1 1) # does not modify image
RegressionTest(TwoD png RD 1) # does not modify image
RegressionTest(FourD nrrd RD)
