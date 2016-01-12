itk_module_test()

set( MorphologicalContourInterpolationTests
  itkMorphologicalContourInterpolationTest.cxx
  )

CreateTestDriver( MorphologicalContourInterpolation "${MorphologicalContourInterpolation-Test_LIBRARIES}" "${MorphologicalContourInterpolationTests}" )

function(JustInputTest ImageName Ext)
itk_add_test( NAME itkMorphologicalContourInterpolationTest_${ImageName}
  COMMAND MorphologicalContourInterpolationTestDriver
  itkMorphologicalContourInterpolationTest
    DATA{Input/${ImageName}.${Ext}}
    ${ITK_TEST_OUTPUT_DIR}/${ImageName}_interpolated.${Ext}
    )
endfunction()

# Handcrafted tests
JustInputTest(SimplestOneToOne nrrd)
JustInputTest(OneToOne nrrd)
JustInputTest(OneToThree nrrd)
JustInputTest(ExtrapolationAppearing nrrd)
JustInputTest(DoubleTwoLabelBranching nrrd)
JustInputTest(TwoAxisDoubleTwoLabelBranchingConflict nrrd)
JustInputTest(ThreeAxisFourLabelConflict nrrd)
JustInputTest(AccidentalMiddleSliceSegmentation nrrd)
JustInputTest(FullEndSevenLabels nrrd)
JustInputTest(ManyToMany nrrd)
JustInputTest(Micro1 nrrd)

# Paul's tests
JustInputTest(105769moving_resliced_to_fixedhz2_segmentation nii.gz)
JustInputTest(105769fixed_segmentation nii.gz)
JustInputTest(64816L_amygdala_fin nii.gz)
JustInputTest(64816L_amygdala_int nii.gz)
JustInputTest(65239R_amygdala_fin nii.gz)
JustInputTest(65239R_amygdala_int nii.gz)
