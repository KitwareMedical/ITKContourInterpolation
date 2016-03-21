/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMorphologicalContourInterpolator.h"
#include <cstdlib>
#include <string>

long int string2int(char *number)
{
  long res = strtol(number, NULL, 10);
  return res;
}

template <typename ImageType>
void doTest(std::string inFilename, std::string outFilename,
  bool UseDistanceTransform, int axis, int label)
{
  typedef itk::ImageFileReader < ImageType > ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(inFilename);

  typedef itk::MorphologicalContourInterpolator<ImageType> mciType;
  typename mciType::Pointer mci = mciType::New();
  mci->SetInput(reader->GetOutput());
  mci->SetUseDistanceTransform(UseDistanceTransform);
  mci->SetAxis(axis);
  mci->SetLabel(label);

  typedef itk::ImageFileWriter< ImageType > WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(outFilename);
  writer->SetInput(mci->GetOutput());
  writer->SetUseCompression(true);
  writer->Update();
}

int itkMorphologicalContourInterpolationTest( int argc, char* argv[] )
{
  if( argc < 3 )
    {
    std::cerr << "Usage: " << argv[0];
    std::cerr << " inputImage outputImage [algo] [axis] [label]\n";
    std::cerr << " algo: RD (repeated dilations ) or DT (distance transform [default])";
    std::cerr << " defaults: axis == -1 (all axes), label == 0 (all labels)";
    std::cerr << std::endl;
    return EXIT_FAILURE;
    }
  const char * inputImageFileName = argv[1];
  const char * outputImageFileName = argv[2];
  bool dt = true; //DistanceTransform
  int axis = -1, label = 0;
  if (argc >= 4)
    {
    std::string algo = argv[3];
    for (auto & c : algo)
      {
      c = toupper(c);
      }
    dt = (algo!="RD");
    }
  if (argc >= 5)
    {
    axis=strtol(argv[4], NULL, 10);
    }
  if (argc >= 6)
    {
    label=strtol(argv[5], NULL, 10);
    }

  typedef itk::ImageIOBase::IOComponentType ScalarPixelType;
  itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(
    inputImageFileName, itk::ImageIOFactory::ReadMode);
  if (!imageIO)
    {
    std::cerr << "Could not CreateImageIO for: " << inputImageFileName << std::endl;
    return EXIT_FAILURE;
    }
  imageIO->SetFileName(inputImageFileName);
  imageIO->ReadImageInformation();
  const ScalarPixelType pixelType = imageIO->GetComponentType();
  const size_t numDimensions = imageIO->GetNumberOfDimensions();

  try
    {
    //unused cases are not instantiated because they greatly increase compile time
    if (numDimensions==2 && pixelType == ScalarPixelType::UCHAR)
      {
      doTest<itk::Image<unsigned char, 2>>(inputImageFileName, outputImageFileName, dt, axis, label);
      return EXIT_SUCCESS;
      }
    if (numDimensions==3 && (pixelType == ScalarPixelType::SHORT
      || pixelType == ScalarPixelType::USHORT))
      {
      doTest<itk::Image<short, 3>>(inputImageFileName, outputImageFileName, dt, axis, label);
      return EXIT_SUCCESS;
      }
    if (numDimensions==4 && pixelType == ScalarPixelType::UCHAR)
      {
      doTest<itk::Image<unsigned char, 4>>(inputImageFileName, outputImageFileName, dt, axis, label);
      return EXIT_SUCCESS;
      }

    std::cerr << "Unsupported image type:\n  Dimensions: " << numDimensions;
    std::cerr<< "\n  Pixel type:" << imageIO->GetComponentTypeAsString(pixelType) << std::endl;
    return EXIT_FAILURE;
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error: " << error << std::endl;
    return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
