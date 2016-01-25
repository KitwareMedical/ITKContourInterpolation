itk_module_test()

set( MorphologicalContourInterpolationTests
  itkMorphologicalContourInterpolationTest.cxx
  )

CreateTestDriver( MorphologicalContourInterpolation "${MorphologicalContourInterpolation-Test_LIBRARIES}" "${MorphologicalContourInterpolationTests}" )

add_executable(runFromIDE manualTest.cxx ${MorphologicalContourInterpolationTests})
target_link_libraries(runFromIDE ${MorphologicalContourInterpolation-Test_LIBRARIES})

function(JustInputTest ImageName Ext)
itk_add_test( NAME itkMorphologicalContourInterpolationTest_${ImageName}
  COMMAND MorphologicalContourInterpolationTestDriver
  itkMorphologicalContourInterpolationTest
    DATA{Input/${ImageName}.${Ext}}
    ${ITK_TEST_OUTPUT_DIR}/${ImageName}_interpolated.${Ext}
    )
endfunction()

# Handcrafted tests
JustInputTest(Micro1 nrrd)
JustInputTest(SimplestOneToOne nrrd)
JustInputTest(OneToOne nrrd)
JustInputTest(OneToThree nrrd)
JustInputTest(ExtrapolationAppearing nrrd)
JustInputTest(DoubleTwoLabelBranching nrrd)
JustInputTest(TwoAxisDoubleTwoLabelBranchingConflict nrrd)
JustInputTest(ThreeAxisFourLabelConflict nrrd) #TODO: debug
JustInputTest(AccidentalMiddleSliceSegmentation nrrd)
JustInputTest(FullEndSevenLabels nrrd) #TODO: debug
JustInputTest(ManyToMany16 nrrd) #TODO: debug
JustInputTest(ManyToMany nrrd) #TODO: debug

# Paul's tests
JustInputTest(105769moving_resliced_to_fixedhz2_segmentation nii.gz)
JustInputTest(105769fixed_segmentation nii.gz)
JustInputTest(64816L_amygdala_int nii.gz)
JustInputTest(64816L_amygdala_fin nii.gz)
JustInputTest(65239R_amygdala_int nii.gz)
JustInputTest(65239R_amygdala_fin nii.gz)
