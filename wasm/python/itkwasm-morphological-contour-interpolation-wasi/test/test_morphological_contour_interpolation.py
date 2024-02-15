from itkwasm_image_io import read_image, write_image

from itkwasm_morphological_contour_interpolation_wasi import morphological_contour_interpolation

from .common import test_input_path, test_output_path

def test_morphological_contour_interpolation():
    test_input_file_path = test_input_path / '64816L_amygdala_int.nii.gz'
    test_output_file_path = test_output_path / '64816L_amygdala_int.nii.gz'

    image = read_image(test_input_file_path)
    interpolated = morphological_contour_interpolation(image)
    write_image(interpolated, test_output_file_path)

    test_output_file_path = test_output_path / '64816L_amygdala_int_axis_1.nii.gz'
    interpolated = morphological_contour_interpolation(image, axis=1)
    write_image(interpolated, test_output_file_path)