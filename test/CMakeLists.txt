itk_module_test()

set( MorphologicalContourInterpolationTests
  itkMorphologicalContourInterpolationTest.cxx
  )

CreateTestDriver( MorphologicalContourInterpolation "${MorphologicalContourInterpolation-Test_LIBRARIES}" "${MorphologicalContourInterpolationTests}" )

add_executable(runFromIDE manualTest.cxx ${MorphologicalContourInterpolationTests})
target_link_libraries(runFromIDE ${MorphologicalContourInterpolation-Test_LIBRARIES})

function(JustInputTest ImageName Ext) #optional: Axis Label
  foreach(f ${ARGN}) #create unique test name
    set(optArgs "${optArgs}_${f}")
  endforeach()
  set(outImage "${ITK_TEST_OUTPUT_DIR}/${ImageName}_interpolated${optArgs}.${Ext}")
  # set(imageNameOpt "${ImageName}${optArgs}") # only used for regression
  
  itk_add_test( NAME itkMorphologicalContourInterpolation_RegressionTest_${ImageName}${optArgs}
    COMMAND MorphologicalContourInterpolationTestDriver
    # --compare DATA{Baseline/${imageNameOpt}.${Ext}} ${outImage} # no regression
    itkMorphologicalContourInterpolationTest
      DATA{Input/${ImageName}.${Ext}} ${outImage} ${ARGN})
endfunction()

function(RegressionTest ImageName Ext) #optional: Axis Label
  foreach(f ${ARGN}) #create unique test name
    set(optArgs "${optArgs}_${f}")
  endforeach()
  set(outImage "${ITK_TEST_OUTPUT_DIR}/${ImageName}_interpolated${optArgs}.${Ext}")
  set(imageNameOpt "${ImageName}${optArgs}")
  
  itk_add_test( NAME itkMorphologicalContourInterpolation_RegressionTest_${ImageName}${optArgs}
    COMMAND MorphologicalContourInterpolationTestDriver
    --compare DATA{Baseline/${imageNameOpt}.${Ext}} ${outImage}
    itkMorphologicalContourInterpolationTest
      DATA{Input/${ImageName}.${Ext}} ${outImage} ${ARGN})
endfunction()

# Handcrafted tests
RegressionTest(Empty nrrd)
RegressionTest(NoSlices nrrd) #interpolation should not be invoked (output==input)
RegressionTest(Micro1 nrrd)
RegressionTest(SimplestOneToOne nrrd)
RegressionTest(OneToOne nrrd)
RegressionTest(OneToThree nrrd)
RegressionTest(ExtrapolationAppearing nrrd)
JustInputTest(DoubleTwoLabelBranching nrrd)
JustInputTest(TwoAxisDoubleTwoLabelBranchingConflict nrrd)
JustInputTest(AccidentalMiddleSliceSegmentation nrrd)
JustInputTest(SevenLabels nrrd)
JustInputTest(FullEnd nrrd)
JustInputTest(ManyToMany16 nrrd)
RegressionTest(ManyToMany nrrd)
RegressionTest(GridSeg nrrd)
RegressionTest(GridSeg2 nrrd)

# Per axis/label tests
RegressionTest(GridSeg nrrd 0)
RegressionTest(GridSeg nrrd 1)
RegressionTest(GridSeg nrrd 2)
RegressionTest(GridSeg2 nrrd 0)
RegressionTest(GridSeg2 nrrd 1)
RegressionTest(GridSeg2 nrrd 2)
RegressionTest(ThreeAxisFourLabelConflict nrrd)
RegressionTest(ThreeAxisFourLabelConflict nrrd 0)
RegressionTest(ThreeAxisFourLabelConflict nrrd 1)
RegressionTest(ThreeAxisFourLabelConflict nrrd 2)
RegressionTest(ThreeAxisFourLabelConflict nrrd 0 0) # equals to _0
RegressionTest(ThreeAxisFourLabelConflict nrrd 1 1)
RegressionTest(ThreeAxisFourLabelConflict nrrd 2 3) # equals to _2
RegressionTest(ThreeAxisFourLabelConflict nrrd 0 4) # equals to _0
RegressionTest(ThreeAxisFourLabelConflict nrrd 1 2)
RegressionTest(ThreeAxisFourLabelConflict nrrd 2 2) # does not modify image


# Paul's tests
JustInputTest(105769moving_resliced_to_fixedhz2_segmentation nii.gz)
JustInputTest(105769fixed_segmentation nii.gz)
JustInputTest(64816L_amygdala_int nii.gz)
JustInputTest(64816L_amygdala_fin nii.gz)
JustInputTest(65239R_amygdala_int nii.gz)
JustInputTest(65239R_amygdala_fin nii.gz)
#TODO: make different combinations of axis and label parameters

# Long-running handcrafted tests
RegressionTest(FaceRoI nrrd)
RegressionTest(FaceRoI nrrd 0)
RegressionTest(FaceRoI nrrd 1)
RegressionTest(FaceRoI nrrd 2)
JustInputTest(BigZ nrrd)
