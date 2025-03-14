#include <iostream>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkMedianImageFilter.h>
#include <itkMorphologicalContourInterpolator.h>


int
main(int argc, char * argv[])
{
  if (argc < 3)
  {
    std::cout << "Syntax: " << argv[0] << " sliceImage interpolatedImage [smoothingRadius=2]" << std::endl;
    return EXIT_FAILURE;
  }

  try
  {
    using ImageType = itk::Image<unsigned char, 3>;
    using ReaderType = itk::ImageFileReader<ImageType>;
    ReaderType::Pointer imageReader = ReaderType::New();
    imageReader->SetFileName(argv[1]);

    using mciType = itk::MorphologicalContourInterpolator<ImageType>;
    mciType::Pointer mci = mciType::New();
    mci->SetInput(imageReader->GetOutput());
    mci->Update();

    int smoothingRadius = 2;
    if (argc >= 4)
    {
      smoothingRadius = std::stoi(argv[3]);
    }

    using MedianType = itk::MedianImageFilter<ImageType, ImageType>;
    MedianType::Pointer medF = MedianType::New();
    medF->SetInput(mci->GetOutput());
    medF->SetRadius(smoothingRadius);
    medF->Update();

    using WriterType = itk::ImageFileWriter<ImageType>;
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
