/*=========================================================================

 *  Copyright NumFOCUS
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         https://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#include "itkInputImage.h"
#include "itkOutputImage.h"
#include "itkPipeline.h"
#include "itkSupportInputImageTypes.h"

#include "itkMorphologicalContourInterpolator.h"

template< typename TImage >
int
MorphologicalContourInterpolation( itk::wasm::Pipeline& pipeline, const TImage* inputImage )
{
  using ImageType = TImage;

  pipeline.get_option( "input-image" )->required()->type_name( "INPUT_IMAGE" );

  using OutputImageType = itk::wasm::OutputImage< ImageType >;
  OutputImageType outputImage;
  pipeline.add_option( "output-image", outputImage, "The output image" )->required()->type_name( "OUTPUT_IMAGE" );

  typename ImageType::IndexValueType label = 0;
  pipeline.add_option( "-l,--label", label, "The label to interpolate. Interpolates all labels if set to 0 (default)." );

  int axis = -1;
  pipeline.add_option( "-a,--axis", axis, "Interpolate only along this axis. Interpolates along all axes if set to -1 (default)." );

  bool noHeuristicAlignment{ false };
  pipeline.add_flag( "--no-heuristic-alignment", noHeuristicAlignment, "Heuristic alignment of regions for interpolation is faster than optimal alignment." );

  bool noUseDistanceTransform{ false };
  pipeline.add_flag( "--no-use-distance-transform",
                     noUseDistanceTransform,
                     "Using distance transform instead of repeated dilations to calculate the median contour is slightly faster, but produces lower quality interpolations." );

  bool useCustomSlicePositions{ false };
  pipeline.add_flag( "--use-custom-slice-positions", useCustomSlicePositions, "Use custom slice positions (not slice auto-detection)." );

  bool noUseExtrapolation{ false };
  pipeline.add_flag( "--no-use-extrapolation",
                     noUseExtrapolation,
                     "Perform extrapolation for branch extremities. Branch extremities are defined as regions having no overlap with any region in the next slice. Extrapolation helps generate smooth "
                     "surface closings." );

  bool useBallStructuringElement;
  pipeline.add_flag( "--use-ball-structuring-element", useBallStructuringElement, "Use ball instead of default cross structuring element for repeated dilations." );

  int labeledSliceIndicesAxis = -1;
  pipeline.add_option( "--labeled-slice-indices-axis", labeledSliceIndicesAxis, "Axis along which the labeled slice indices are defined. Default is -1 (that is, auto-detection)." );

  int labeledSliceIndicesLabel = 1;
  pipeline.add_option( "--labeled-slice-indices-label", labeledSliceIndicesLabel, "Label of the slice indices. Default is 1." );

  std::vector< typename ImageType::IndexValueType > labeledSliceIndices;
  pipeline.add_option( "--labeled-slice-indices", labeledSliceIndices, "List of labeled slice indices. Default is empty." );

  ITK_WASM_PARSE( pipeline );

  using InterpolatorType = itk::MorphologicalContourInterpolator< ImageType >;
  auto interpolator = InterpolatorType::New();
  interpolator->SetInput( inputImage );

  interpolator->SetLabel( label );
  interpolator->SetAxis( axis );
  interpolator->SetHeuristicAlignment( !noHeuristicAlignment );
  interpolator->SetUseDistanceTransform( !noUseDistanceTransform );
  interpolator->SetUseCustomSlicePositions( useCustomSlicePositions );
  interpolator->SetUseExtrapolation( !noUseExtrapolation );
  interpolator->SetUseBallStructuringElement( !useBallStructuringElement );
  if ( labeledSliceIndicesAxis != -1 && !labeledSliceIndices.empty() )
    {
      interpolator->SetLabeledSliceIndices( labeledSliceIndicesAxis, labeledSliceIndicesLabel, labeledSliceIndices );
    }

  ITK_WASM_CATCH_EXCEPTION( pipeline, interpolator->UpdateLargestPossibleRegion() );

  typename ImageType::ConstPointer output = interpolator->GetOutput();
  outputImage.Set( output );

  return EXIT_SUCCESS;
}

template< typename TImage >
class PipelineFunctor
{
public:
  int
  operator()( itk::wasm::Pipeline& pipeline )
  {
    using ImageType = TImage;

    using InputImageType = itk::wasm::InputImage< ImageType >;
    InputImageType inputImage;
    pipeline.add_option( "input-image", inputImage, "The input image" );

    ITK_WASM_PRE_PARSE( pipeline );

    typename ImageType::ConstPointer image = inputImage.Get();
    return MorphologicalContourInterpolation< ImageType >( pipeline, image );
  }
};

int
main( int argc, char* argv[] )
{
  itk::wasm::Pipeline pipeline( "morphological-contour-interpolation", "Interpolates contours between slices.", argc, argv );

  return itk::wasm::SupportInputImageTypes< PipelineFunctor, uint8_t, int8_t, uint16_t, int16_t, uint32_t, int32_t, uint64_t, int64_t >::Dimensions< 3U, 4U >( "input-image", pipeline );
}
