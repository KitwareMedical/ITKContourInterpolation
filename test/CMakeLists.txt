itk_module_test()

set( MorphologicalContourInterpolationTests
  itkMorphologicalContourInterpolationTest.cxx
  )

CreateTestDriver( MorphologicalContourInterpolation "${MorphologicalContourInterpolation-Test_LIBRARIES}" "${MorphologicalContourInterpolationTests}" )

add_executable(runFromIDE manualTest.cxx ${MorphologicalContourInterpolationTests})
target_link_libraries(runFromIDE ${MorphologicalContourInterpolation-Test_LIBRARIES})

function(JustInputTest ImageName Ext)
itk_add_test( NAME itkMorphologicalContourInterpolation_JustInputTest_${ImageName}
  COMMAND MorphologicalContourInterpolationTestDriver
  itkMorphologicalContourInterpolationTest
    DATA{Input/${ImageName}.${Ext}}
    ${ITK_TEST_OUTPUT_DIR}/${ImageName}_interpolated.${Ext}
    )
endfunction()

function(RegressionTest ImageName Ext)
itk_add_test( NAME itkMorphologicalContourInterpolation_RegressionTest_${ImageName}
  COMMAND MorphologicalContourInterpolationTestDriver
  --compare DATA{Baseline/${ImageName}.${Ext}}
            ${ITK_TEST_OUTPUT_DIR}/${ImageName}_interpolated.${Ext}
  itkMorphologicalContourInterpolationTest
    DATA{Input/${ImageName}.${Ext}}
    ${ITK_TEST_OUTPUT_DIR}/${ImageName}_interpolated.${Ext}
    )
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
JustInputTest(ThreeAxisFourLabelConflict nrrd)
JustInputTest(AccidentalMiddleSliceSegmentation nrrd)
JustInputTest(SevenLabels nrrd)
JustInputTest(FullEnd nrrd)
JustInputTest(ManyToMany16 nrrd)
RegressionTest(ManyToMany nrrd)
JustInputTest(GridSeg nrrd)

# Paul's tests
JustInputTest(105769moving_resliced_to_fixedhz2_segmentation nii.gz)
JustInputTest(105769fixed_segmentation nii.gz)
JustInputTest(64816L_amygdala_int nii.gz)
JustInputTest(64816L_amygdala_fin nii.gz)
JustInputTest(65239R_amygdala_int nii.gz)
JustInputTest(65239R_amygdala_fin nii.gz)

# Long-running handcrafted tests
JustInputTest(BigZ nrrd)
JustInputTest(FaceRoI nrrd) #TODO: convert into regression test with FaceRoIFinished.nrrd as expected output
