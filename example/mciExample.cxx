#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkMedianImageFilter.h>
#include <itkMorphologicalContourInterpolator.h>
#include <iostream>


int main(int argc, char * argv[])
{
  if (argc < 3)
    {
    std::cout << "Syntax: " << argv[0] << " sliceImage interpolatedImage [smoothingRadius=2]" << std::endl;
    return EXIT_FAILURE;
    }

  try
    {
    typedef itk::Image<unsigned char, 3> ImageType;
    typedef itk::ImageFileReader<ImageType> ReaderType;
    ReaderType::Pointer imageReader = ReaderType::New();
    imageReader->SetFileName(argv[1]);

    typedef itk::MorphologicalContourInterpolator<ImageType> mciType;
    mciType::Pointer mci = mciType::New();
    mci->SetInput(imageReader->GetOutput());
    mci->Update();

    int smoothingRadius = 2;
    if (argc >= 4)
      {
      smoothingRadius = atoi(argv[3]);
      }

    typedef itk::MedianImageFilter<ImageType, ImageType> MedianType;
    MedianType::Pointer medF = MedianType::New();
    medF->SetInput(mci->GetOutput());
    medF->SetRadius(smoothingRadius);
    medF->Update();

    typedef itk::ImageFileWriter<ImageType> WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(argv[2]);
    writer->SetInput(medF->GetOutput());
    writer->SetUseCompression(true);
    writer->Update();
    }
  catch (itk::ExceptionObject & exc)
    {
    std::cerr << exc << std::endl;
    return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}