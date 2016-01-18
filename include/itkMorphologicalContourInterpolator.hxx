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
#ifndef itkMorphologicalContourInterpolator_hxx
#define itkMorphologicalContourInterpolator_hxx

#include "itkMorphologicalContourInterpolator.h"
#include "itkObjectFactory.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageAlgorithm.h"
#include "itkUnaryFunctorImageFilter.h"
#include "itkCastImageFilter.h"
#include <utility>
#include <algorithm>
#include <queue>

//DEBUG
#include <iostream>
#include "itkImageFileWriter.h"


namespace itk
{
template< typename TImage >
void WriteDebug(typename TImage::Pointer out, const char *filename)
{
  typedef ImageFileWriter<TImage> WriterType;
  typename WriterType::Pointer w = WriterType::New();
  w->SetInput(out);
  w->SetFileName(filename);
  try
    {
    w->Update();
    }
  catch (itk::ExceptionObject & error)
    {
    std::cerr << "Error: " << error << std::endl;
    }
}

void WriteDebug(itk::Image<bool, 3>::Pointer out, const char *filename)
{
  typedef itk::Image<bool, 3>                                 BoolImageType;
  typedef itk::Image<unsigned char, 3>                        ucharImageType;
  typedef itk::CastImageFilter<BoolImageType, ucharImageType> CastType;
  CastType::Pointer caster = CastType::New();
  caster->SetInput(out);
  WriteDebug<ucharImageType>(caster->GetOutput(), filename);
}

template< typename TImage >
bool ImagesEqual(typename TImage::Pointer a, typename TImage::Pointer b)
{
  ImageRegionConstIterator<TImage> ita(a, a->GetLargestPossibleRegion());
  ImageRegionConstIterator<TImage> itb(b, b->GetLargestPossibleRegion());

  while (!ita.IsAtEnd())
    {
    if (ita.Get() != itb.Get())
      {
      break;
      }
    ++ita;
    ++itb;
    }

  if (ita.IsAtEnd())
    {
    return true;
    }
  else
    {
    return false;
    }
}


template< typename TImage >
MorphologicalContourInterpolator<TImage>
::MorphologicalContourInterpolator()
 :m_Label(0),
  m_Axis(-1),
  m_HeuristicAlignment(true),
  m_LabeledSlices(TImage::ImageDimension) //initialize with empty sets
{
  m_Or = OrType::New();
  m_Dilator = DilateType::New(); //structuring element in InterpolateAlong
  m_And = AndFilterType::New();

  //set up pipeline for regioned connected components
  m_RoI = RoiType::New();
  m_Binarizer = BinarizerType::New();
  m_Binarizer->SetInput(m_RoI->GetOutput());
  m_ConnectedComponents = ConnectedComponentsType::New();
  m_ConnectedComponents->SetInput(m_Binarizer->GetOutput());
  m_ConnectedComponents->SetFullyConnected(true);
}


template< typename TImage >
typename MorphologicalContourInterpolator<TImage>::SliceSetType
MorphologicalContourInterpolator<TImage>
::GetLabeledSliceIndices(unsigned int axis)
{
  return m_LabeledSlices[axis];
}


template< typename TImage >
void
MorphologicalContourInterpolator<TImage>
::SetLabeledSliceIndices(unsigned int axis, SliceSetType indices)
{
  m_LabeledSlices[axis]=indices;
  this->Modified();
}


template< typename TImage >
void
MorphologicalContourInterpolator<TImage>
::SetLabeledSliceIndices(unsigned int axis, std::vector<typename TImage::IndexValueType> indices)
{
  m_LabeledSlices[axis] = SliceSetType().insert(indices.begin(), indices.end());
  this->Modified();
}


template< typename TImage >
void
MorphologicalContourInterpolator<TImage>
::ExpandRegion(typename TImage::RegionType &region, typename TImage::IndexType index)
{
  for (unsigned int a = 0; a < TImage::ImageDimension; ++a)
    {
    if (region.GetIndex(a) > index[a])
      {
      region.SetSize(a, region.GetSize(a) + region.GetIndex(a) - index[a]);
      region.SetIndex(a, index[a]);
      }
    else if (region.GetIndex(a) + region.GetSize(a) <= index[a])
      {
      region.SetSize(a, index[a] - region.GetIndex(a) + 1);
      }
    //else it is already within
    }
}


template< typename TImage >
void
MorphologicalContourInterpolator<TImage>
::DetermineSliceOrientations()
{
  m_LabeledSlices.clear();
  m_LabeledSlices.resize(TImage::ImageDimension); //initialize with empty sets
  m_BoundingBoxes.clear();
  m_Orientations.clear();

  typename TImage::RegionType region = m_Output->GetRequestedRegion();
  ImageRegionConstIteratorWithIndex<TImage> it(m_Input, region);

  OrientationType orientations = OrientationType();
  orientations.Fill(false);

  for (; !it.IsAtEnd(); ++it)
    {
    typename TImage::IndexType indPrev, indNext;
    const typename TImage::IndexType ind = it.GetIndex();
    const typename TImage::PixelType val = m_Input->GetPixel(ind);
    if (val != 0 || (m_Label != 0 && val == m_Label))
      {
      typename TImage::RegionType boundingBox1;
      boundingBox1.SetIndex(ind);
      for (unsigned int a = 0; a < TImage::ImageDimension; ++a)
        {
        boundingBox1.SetSize(a, 1);
        }
      std::pair<typename BoundingBoxesType::iterator, bool> resBB
        = m_BoundingBoxes.insert(std::make_pair(val, boundingBox1));
      if (!resBB.second) //include this index in existing BB
        {
        ExpandRegion(resBB.first->second, ind);
        }

      std::pair<typename OrientationsType::iterator, bool> res
        = m_Orientations.insert(std::make_pair(val, orientations));
      typename OrientationsType::iterator oRef = res.first;
      unsigned int cTrue = 0;
      unsigned int cAdjacent = 0;
      unsigned int axis = 0;
      for (unsigned int a = 0; a < TImage::ImageDimension; ++a)
        {
        indPrev = ind;
        indPrev[a]--;
        indNext = ind;
        indNext[a]++;
        typename TImage::PixelType prev = 0;
        if (region.IsInside(indPrev))
          {
          prev = m_Input->GetPixel(indPrev);
          }
        typename TImage::PixelType next = 0;
        if (region.IsInside(indNext))
          {
          next = m_Input->GetPixel(indNext);
          }
        if (prev == 0 && next == 0)
          {
          axis = a;
          ++cTrue;
          }
        else if (prev == val && next == val)
          {
          ++cAdjacent;
          }
        }
      if (cTrue == 1 && cAdjacent == TImage::ImageDimension - 1) //slice has empty adjacent space only along one axis
        {
        oRef->second[axis] = true; //add this dimension for this label
        if (m_Axis == -1 || m_Axis == axis)
          {
          m_LabeledSlices[axis][val].insert(ind[axis]);
          }
        }
      }
    }
}


template< typename TImage >
void
MorphologicalContourInterpolator<TImage>
::Extrapolate(int axis, TImage *out, typename TImage::PixelType label,
  typename TImage::IndexValueType i, typename TImage::IndexValueType j,
  typename TImage::Pointer iConn, typename TImage::PixelType iRegionId)
{
  PixelList jRegionIds;
  jRegionIds.push_back(iRegionId);
  typename TImage::IndexType centroid = Centroid(iConn, jRegionIds);
  centroid[axis] = j;

  typename TImage::RegionType reg3;
  typename TImage::SizeType size;
  size.Fill(3);
  size[axis] = 1;
  reg3.SetSize(size);

  typename TImage::IndexType phIndex;
  for (unsigned d = 0; d < TImage::ImageDimension; d++)
    {
    phIndex[d] = centroid.GetIndex()[d] - 1;
    }
  phIndex[axis] = j;
  reg3.SetIndex(phIndex);

  //create a phantom small slice centered around centroid
  typename TImage::Pointer phSlice = TImage::New();
  phSlice->CopyInformation(iConn);
  phSlice->SetRegions(reg3);
  phSlice->Allocate(true);

  //add a phantom point to the center of a newly constructed slice
  phSlice->SetPixel(centroid, 1);
  //WriteDebug<TImage>(iConn, "C:\\iConn.nrrd");
  //WriteDebug<TImage>(phSlice, "C:\\phSlice.nrrd");
  typename TImage::IndexType translation = Align(axis, iConn, iRegionId, phSlice, jRegionIds);
  Interpolate1to1(axis, out, label, i, j, iConn, iRegionId, phSlice, 1, translation);
}


template< typename TImage >
typename MorphologicalContourInterpolator<TImage>::BoolImageType::Pointer
MorphologicalContourInterpolator<TImage>
::Dilate1(typename BoolImageType::Pointer seed, typename BoolImageType::Pointer mask)
{
  m_Dilator->SetInput(seed);
  m_Dilator->GetOutput()->SetRegions(seed->GetRequestedRegion());
  m_Dilator->Update();
  typename BoolImageType::Pointer temp = m_Dilator->GetOutput();
  temp->DisconnectPipeline();
  //temp->SetRegions(mask->GetLargestPossibleRegion()); //not needed when seed and mask have same regions
  m_And->SetInput(0, mask);
  m_And->SetInput(1, temp);
  m_And->GetOutput()->SetRegions(seed->GetRequestedRegion());
  m_And->Update();
  typename BoolImageType::Pointer result = m_And->GetOutput();
  result->DisconnectPipeline();
  WriteDebug(seed, "C:\\seed.nrrd");
  WriteDebug(mask, "C:\\mask.nrrd");
  WriteDebug(temp, "C:\\temp.nrrd");
  WriteDebug(result, "C:\\result.nrrd");
  return result;
}


template< typename TImage >
std::vector<typename MorphologicalContourInterpolator<TImage>::BoolImageType::Pointer>
MorphologicalContourInterpolator<TImage>
::GenerateDilationSequence(
  typename BoolImageType::Pointer begin, typename BoolImageType::Pointer end)
{
  std::vector<typename BoolImageType::Pointer> seq;
  seq.push_back(Dilate1(begin, end));
  do
    {
    seq.back()->DisconnectPipeline();
    seq.push_back(Dilate1(seq.back(), end));
    } while (!ImagesEqual<BoolImageType>(seq.back(), seq[seq.size() - 2]));
  seq.pop_back(); //remove duplicate image
  return seq;
}


template< typename TImage >
void
MorphologicalContourInterpolator<TImage>
::Interpolate1to1(int axis, TImage *out, typename TImage::PixelType label,
  typename TImage::IndexValueType i, typename TImage::IndexValueType j,
  typename TImage::Pointer iConn, typename TImage::PixelType iRegionId,
  typename TImage::Pointer jConn, typename TImage::PixelType jRegionId,
  typename TImage::IndexType translation)
{
  //translate iConn by t/2 and jConn by -t/2
  typename TImage::IndexType iTrans;
  typename TImage::IndexType jTrans;
  typename TImage::RegionType iRegion = iConn->GetLargestPossibleRegion();
  typename TImage::RegionType jRegion = jConn->GetLargestPossibleRegion();
  typename TImage::RegionType newRegion;
  //translation[axis] = j - i;
  for (unsigned d = 0; d < TImage::ImageDimension; d++)
    {
    iTrans[d] = translation[d] / 2;
    jTrans[d] = iTrans[d] - translation[d];
    newRegion.SetIndex(d, (iRegion.GetIndex()[d] + jRegion.GetIndex()[d]) / 2);
    newRegion.SetSize(d, std::max(iRegion.GetSize(d), jRegion.GetSize(d)));
    }
  //newRegion.SetIndex(axis, (i + j) / 2); //redundant
  typename TImage::Pointer iConnT = TranslateImage(iConn, iTrans, newRegion);
  typename TImage::Pointer jConnT = TranslateImage(jConn, jTrans, newRegion);
  WriteDebug<TImage>(iConnT, "C:\\iConnT.nrrd");
  WriteDebug<TImage>(jConnT, "C:\\jConnT.nrrd");

  //convert to binary masks
  MatchesID matchesIDi(iRegionId);
  MatchesID matchesIDj(jRegionId);
  typedef itk::UnaryFunctorImageFilter<TImage, BoolImageType, MatchesID> CastType;
  typename CastType::Pointer caster = CastType::New();
  caster->SetFunctor(matchesIDi);
  caster->SetInput(iConnT);
  caster->Update();
  typename BoolImageType::Pointer iMask = caster->GetOutput();
  iMask->DisconnectPipeline();
  caster->SetFunctor(matchesIDj);
  caster->SetInput(jConnT);
  caster->Update();
  typename BoolImageType::Pointer jMask = caster->GetOutput();
  jMask->DisconnectPipeline();

  //generate sequence
  WriteDebug<TImage>(iConn, "C:\\iConn.nrrd");
  WriteDebug<TImage>(jConn, "C:\\jConn.nrrd");
  WriteDebug(iMask, "C:\\iMask.nrrd");
  WriteDebug(jMask, "C:\\jMask.nrrd");
  m_And->SetInput(0, iMask);
  m_And->SetInput(1, jMask);
  m_And->GetOutput()->SetRegions(newRegion);
  m_And->Update();
  typename BoolImageType::Pointer intersection = m_And->GetOutput();
  intersection->DisconnectPipeline();
  WriteDebug(intersection, "C:\\intersection.nrrd");
  std::vector<typename BoolImageType::Pointer> iSeq = GenerateDilationSequence(intersection, iMask);
  std::vector<typename BoolImageType::Pointer> jSeq = GenerateDilationSequence(intersection, jMask);
  unsigned minSeq = std::min(iSeq.size(), jSeq.size());
  unsigned maxSeq = std::max(iSeq.size(), jSeq.size());
  std::reverse(iSeq.begin(), iSeq.end()); //we want to start from i and end at intersection
  std::vector<typename BoolImageType::Pointer> seq; 
  for (unsigned x = 0; x < minSeq; x++)
    {
    m_Or->SetInput(0, iSeq[x]);
    m_Or->SetInput(1, jSeq[x]);
    m_Or->GetOutput()->SetRegions(newRegion);
    m_Or->Update();
    seq.push_back(m_Or->GetOutput());
    seq.back()->DisconnectPipeline();
    }
  if (iSeq.size() < jSeq.size())
    {
    for (unsigned x = minSeq; x < maxSeq; x++)
      {
      seq.push_back(jSeq[x]);
      }
    }
  else
    {
    for (unsigned x = minSeq; x < maxSeq; x++)
      {
      seq.push_back(iSeq[x]);
      }
    }

  //find median
  unsigned minIndex;
  IdentifierType min = newRegion.GetNumberOfPixels();
  for (unsigned x = 0; x < maxSeq; x++)
    {
    WriteDebug(seq[x], (std::string("C:\\seq") + char('A' + x) + ".nrrd").c_str());
    IdentifierType iS = CardSymDifference(seq[x], iMask);
    IdentifierType jS = CardSymDifference(seq[x], jMask);
    IdentifierType xScore = iS >= jS ? iS - jS : jS - iS; //abs(iS-jS)
    if (xScore < min)
      {
      min = xScore;
      minIndex = x;
      }
    }

  //finally write it out into the output image pointer
  typename TImage::RegionType outRegion = m_Input->GetLargestPossibleRegion();
  typename TImage::IndexType t0 = { 0 };
  IntersectionRegions(t0, outRegion,newRegion);
  ImageRegionConstIterator<BoolImageType> seqIt(seq[minIndex], newRegion);
  ImageRegionIterator<TImage> outIt(out, newRegion);
  while (!outIt.IsAtEnd())
    {
    if (seqIt.Get())
      {
      outIt.Set(label);
      }
    ++outIt;
    ++seqIt;
    }

  //recurse if needed
  typename TImage::IndexValueType mid = (i + j) / 2;
  if (abs(i - j) > 2)
    {
      typename TImage::Pointer midConn = TImage::New();
      midConn->CopyInformation(seq[minIndex]);
      midConn->SetRegions(newRegion);
      midConn->Allocate(true);
      ImageAlgorithm::Copy<TImage, TImage>(out, midConn.GetPointer(), newRegion, newRegion);
      WriteDebug<TImage>(out, "C:\\midConnSource.nrrd");
      WriteDebug<TImage>(midConn, "C:\\midConn.nrrd");
      //seqIt.GoToBegin();
      //ImageRegionIterator<TImage> outIt(midConn, newRegion);
      //while (!outIt.IsAtEnd())
      //  {
      //  if (seqIt.Get())
      //    {
      //    outIt.Set(label);
      //    }
      //  ++outIt;
      //  ++seqIt;
      //  }
      PixelList regionIDs;
      regionIDs.push_back(label);

      if (abs(i - mid) > 1)
        {
        typename TImage::IndexType tRecurse = Align(axis, iConn, iRegionId, midConn, regionIDs);
        //optimization tRecurse=iTrans
        Interpolate1to1(axis, out, label, i, mid, iConn, iRegionId, midConn, label, tRecurse);
        }
      if (abs(j - mid) > 1)
        {
        typename TImage::IndexType tRecurse = Align(axis, jConn, jRegionId, midConn, regionIDs);
        //optimization tRecurse=jTrans
        Interpolate1to1(axis, out, label, j, mid, jConn, jRegionId, midConn, label, tRecurse);
        }
    }
}


template< typename TImage >
void
MorphologicalContourInterpolator<TImage>
::Interpolate1toN(int axis, TImage *out, typename TImage::PixelType label,
typename TImage::IndexValueType i, typename TImage::IndexValueType j,
typename TImage::Pointer iConn, typename TImage::PixelType iRegionId,
typename TImage::Pointer jConn, PixelList jRegionIds,
typename TImage::IndexType translation)
{
  //first convert iConn into binary mask
  MatchesID matchesID(iRegionId);
  
  typedef itk::UnaryFunctorImageFilter<TImage, BoolImageType, MatchesID> CastType;
  typename CastType::Pointer caster = CastType::New();
  caster->SetFunctor(matchesID);
  caster->SetInput(iConn);
  caster->Update();
  typename BoolImageType::Pointer mask = caster->GetOutput();
  //WriteDebug(mask, "C:\\mask.nrrd");
  //WriteDebug<TImage>(iConn, "C:\\iConn.nrrd");
  //WriteDebug<TImage>(jConn, "C:\\jConn.nrrd");

  typename TImage::RegionType iRegion, jRegion, newjRegion;
  iRegion = iConn->GetLargestPossibleRegion();
  jRegion = jConn->GetLargestPossibleRegion();
  newjRegion = jRegion;
  newjRegion.SetSize(iRegion.GetSize());

  //construct n empty images
  std::vector<typename BoolImageType::Pointer> blobs;
  for (unsigned x = 0; x < jRegionIds.size(); x++)
    {
    typename BoolImageType::Pointer temp = BoolImageType::New();
    temp->CopyInformation(jConn);
    temp->SetRegions(iRegion);
    temp->Allocate(true);
    blobs.push_back(temp);
    }

  //fill the n images with intersections - these are seeds
  typename TImage::Pointer belongs = TImage::New();
  belongs->CopyInformation(mask);
  belongs->SetRegions(iRegion);
  belongs->Allocate(true); //initialize to zero (false)
  ImageRegionIterator<TImage> belongIt(belongs, iRegion);
  IntersectionRegions(translation, iRegion, jRegion);
  ImageRegionConstIterator<BoolImageType> maskIt(mask, iRegion);
  ImageRegionConstIteratorWithIndex<TImage> jIt(jConn, jRegion);
  ImageRegionIterator<TImage> belongInit(belongs, iRegion);

  //convert jConn into n blobs, translating them into the index space of iConn
  while (!maskIt.IsAtEnd())
    {
    if (maskIt.Get())
      {
      typename TImage::PixelType jVal=jIt.Get();
      typename PixelList::iterator res = std::find(jRegionIds.begin(), jRegionIds.end(), jVal);
      if (res != jRegionIds.end())
        {
        blobs[res - jRegionIds.begin()]->SetPixel(maskIt.GetIndex(), true);
        belongInit.Set(res - jRegionIds.begin() + 1);
        }
      }
    ++maskIt;
    ++jIt;
    ++belongInit;
    }
  WriteDebug<TImage>(belongs, "C:\\belongs.nrrd");

  //prepare dilation filter
  iRegion = iConn->GetLargestPossibleRegion(); //expand to full i image
  for (unsigned x = 0; x < jRegionIds.size(); x++)
    {
    blobs[x]->SetRegions(iRegion);
    WriteDebug(blobs[x], (std::string("C:\\blob") + char('0' + x) + ".nrrd").c_str());
    }
  ImageRegionConstIterator<BoolImageType> maskIt2(mask, iRegion);
  ImageRegionConstIteratorWithIndex<BoolImageType> jIt2(blobs[0], iRegion);
  
  bool hollowedMaskEmpty;
  do //while hollowed mask is not empty
    {
    for (unsigned x = 0; x < jRegionIds.size(); x++)
      {
      blobs[x] = Dilate1(blobs[x], mask);
      WriteDebug(blobs[x], (std::string("C:\\blob") + char('0' + x) + ".nrrd").c_str());
      blobs[x]->DisconnectPipeline();
      //TODO: save these in a sequence so we don't have to recalculate it!
      }

    hollowedMaskEmpty = true;
    maskIt2.GoToBegin();
    jIt2.GoToBegin();
    belongIt.GoToBegin();
    while (!maskIt2.IsAtEnd()) //hollow out the big mask with dilated seeds while avoiding conflicts
      {
      if (maskIt2.Get())
        {
        if (!belongIt.Get())
          {
          unsigned x = 0;
          for (; x < jRegionIds.size(); x++)
            {
            if (blobs[x]->GetPixel(jIt2.GetIndex()))
              {
              break;
              }
            }
          if (x < jRegionIds.size()) //covered by a blob, hollow it out
            {
            belongIt.Set(x + 1);
            for (x++; x < jRegionIds.size(); x++)
              {
              //pixel does not belong to this blob
              blobs[x]->SetPixel(jIt2.GetIndex(), false);
              }
            }
          else //keep it
            {
            hollowedMaskEmpty = false;
            }
          }
        else //the pixel already belongs to some blob
          {
          for (unsigned x = 0; x < jRegionIds.size(); x++)
            {
            if (belongIt.Get() != x + 1)
              {
              //pixel does not belong to this blob
              blobs[x]->SetPixel(jIt2.GetIndex(), false);
              }
            }
          }
        }
      ++maskIt2;
      ++jIt2;
      ++belongIt;
      }
    WriteDebug<TImage>(belongs, "C:\\belongs.nrrd");
    } while (!hollowedMaskEmpty);
  blobs.clear(); // deallocates the images

  //convert the belongs into n Conn-style images
  std::vector<typename TImage::Pointer> conns;
  for (unsigned x = 0; x < jRegionIds.size(); x++)
    {
    typename TImage::Pointer temp2 = TImage::New();
    temp2->CopyInformation(iConn);
    temp2->SetRegions(iConn->GetLargestPossibleRegion());
    temp2->Allocate(true);
    conns.push_back(temp2);
    }
  ImageRegionConstIteratorWithIndex<TImage> belongIt2(belongs, iRegion);
  while (!belongIt2.IsAtEnd())
    {
    const typename TImage::PixelType belong = belongIt2.Get();
    if (belong > 0)
      {
      conns[belong - 1]->SetPixel(belongIt2.GetIndex(), iRegionId);
      }
    ++belongIt2;
    }

  //make n 1-to-1 interpolations
  for (unsigned x = 0; x < jRegionIds.size(); x++)
    {
    typename TImage::IndexType t0 = { 0 };
    typename TImage::Pointer temp = TranslateImage(conns[x], t0, iRegion);
    Interpolate1to1(axis, out, label, i, j, temp, iRegionId, jConn, jRegionIds[x], translation);
    //TODO: call sequence construction directly from here!
    }
}


template< typename TImage >
typename TImage::RegionType
MorphologicalContourInterpolator<TImage>
::MergeBoundingBoxes(const BoundingBoxesType& boundingBoxes)
{
  typename BoundingBoxesType::iterator it = m_BoundingBoxes.begin();
  typename TImage::RegionType result = it->second;
  typename TImage::SizeType minusOne;
  minusOne.Fill(-1);
  for (++it; it != m_BoundingBoxes.end(); ++it)
    {
    ExpandRegion(result, it->second.GetIndex());
    ExpandRegion(result, it->second.GetIndex() + it->second.GetSize() + minusOne);
    }
  return result;
}


template< typename TImage >
typename TImage::Pointer
MorphologicalContourInterpolator<TImage>
::TranslateImage(typename TImage::Pointer image,
  typename TImage::IndexType translation, typename TImage::RegionType newRegion)
{
    //needed?
  typename TImage::Pointer result = TImage::New();
  result->CopyInformation(image);
  result->SetRegions(newRegion);
  result->Allocate(true); //initialize to zero (false)
  typename TImage::RegionType inRegion = image->GetLargestPossibleRegion();
  typename TImage::RegionType rRegion = image->GetLargestPossibleRegion();
  IntersectionRegions(translation, inRegion, newRegion);
  ImageRegionIterator<TImage> resultIt(result, newRegion);
  ImageRegionConstIterator<TImage> inIt(image, inRegion);

  while (!resultIt.IsAtEnd())
    {
    resultIt.Set(inIt.Get());
    ++resultIt;
    ++inIt;
    }
  return result;
}


template< typename TImage >
void MorphologicalContourInterpolator<TImage>
::IntersectionRegions(typename TImage::IndexType translation,
  typename TImage::RegionType & iRegion, typename TImage::RegionType & jRegion)
{
  typename TImage::IndexType iBegin = iRegion.GetIndex();
  typename TImage::IndexType jBegin = jRegion.GetIndex();
  for (IdentifierType d = 0; d < TImage::ImageDimension; d++)
    {
    IdentifierType iSize = iRegion.GetSize(d);
    IdentifierType jSize = jRegion.GetSize(d);
    jBegin[d] -= translation[d];
    IndexValueType t = std::max(iBegin[d], jBegin[d]);
    iRegion.SetSize(d, std::min(iSize - (t - iBegin[d]), jSize - (t - jBegin[d])));
    iRegion.SetIndex(d, t);
    jRegion.SetIndex(d, t + translation[d]);
    }
  jRegion.SetSize(iRegion.GetSize()); //size is the same
}


template< typename TImage >
IdentifierType MorphologicalContourInterpolator<TImage>
::Intersection(typename TImage::Pointer iConn, typename TImage::PixelType iRegionId,
  typename TImage::Pointer jConn, PixelList jRegionIds,
  typename TImage::IndexType translation)
{
  typename TImage::RegionType iRegion, jRegion;
  iRegion = iConn->GetLargestPossibleRegion();
  jRegion = jConn->GetLargestPossibleRegion();
  IntersectionRegions(translation, iRegion, jRegion);

  IdentifierType count = 0;
  ImageRegionConstIterator<TImage> iIt(iConn, iRegion);
  ImageRegionConstIterator<TImage> jIt(jConn, jRegion);
  while (!iIt.IsAtEnd())
    {
    if (iIt.Get() == iRegionId)
      {
      typename TImage::PixelType jVal=jIt.Get();
      typename PixelList::iterator res = std::find(jRegionIds.begin(), jRegionIds.end(), jVal);
      if (res != jRegionIds.end())
        {
        count++;
        }
      }
    ++iIt;
    ++jIt;
    }
  return count;
}

template< typename TImage >
IdentifierType MorphologicalContourInterpolator<TImage>
::CardSymDifference(typename MorphologicalContourInterpolator<TImage>::BoolImageType::Pointer iShape,
  typename MorphologicalContourInterpolator<TImage>::BoolImageType::Pointer jShape)
{
  typename TImage::RegionType region = iShape->GetLargestPossibleRegion();
  IdentifierType count = 0;
  ImageRegionConstIterator<BoolImageType> iIt(iShape, region);
  ImageRegionConstIterator<BoolImageType> jIt(jShape, region);
  while (!iIt.IsAtEnd())
    {
    if (iIt.Get() != jIt.Get())
      {
      count++;
      }
    ++iIt;
    ++jIt;
    }
  return count;
}

template< typename TImage >
typename TImage::IndexType
MorphologicalContourInterpolator<TImage>
::Centroid(typename TImage::Pointer conn, PixelList regionIds)
{
  ImageRegionConstIteratorWithIndex<TImage> it(conn, conn->GetLargestPossibleRegion());
  IndexValueType ind[TImage::ImageDimension] = { 0 }; //all components are initialized to zero
  IdentifierType pixelCount = 0;
  while (!it.IsAtEnd())
    {
    typename TImage::PixelType val = it.Get();
    if (val)
      {
      typename PixelList::iterator res = std::find(regionIds.begin(), regionIds.end(), val);
      if (res != regionIds.end())
        {
        ++pixelCount;
        typename TImage::IndexType pInd = it.GetIndex();
        for (unsigned d = 0; d < TImage::ImageDimension; d++)
          {
          ind[d] += pInd[d];
          }
        }
      }
    ++it;
    }
  typename TImage::IndexType retVal;
  for (unsigned d = 0; d < TImage::ImageDimension; d++)
    {
    retVal[d] = ind[d] / pixelCount;
    }
  return retVal;
}


template< typename TImage >
typename TImage::IndexType
MorphologicalContourInterpolator<TImage>
::Align(int axis, typename TImage::Pointer iConn, typename TImage::PixelType iRegionId,
typename TImage::Pointer jConn, PixelList jRegionIds)
{
  //calculate centroids
  PixelList iRegionIds;
  iRegionIds.push_back(iRegionId);
  typename TImage::IndexType iCentroid = Centroid(iConn, iRegionIds);
  typename TImage::IndexType jCentroid = Centroid(jConn, jRegionIds);

  typename TImage::IndexType ind, centroidInd;
  for (unsigned d = 0; d < TImage::ImageDimension; d++)
    {
    ind[d] = jCentroid[d] - iCentroid[d];
    }
  //ind[axis] = 0; //i and j have different coordinate along this axis
  centroidInd = ind;

  //construct an image with all possible translations
  typename TImage::RegionType searchRegion;
  typename TImage::RegionType iLPR = iConn->GetLargestPossibleRegion();
  typename TImage::RegionType jLPR = jConn->GetLargestPossibleRegion();
  for (IdentifierType d = 0; d < TImage::ImageDimension; d++)
    {
    searchRegion.SetIndex(d, /*jLPR.GetIndex()[d]*/ - iLPR.GetIndex()[d] - jLPR.GetSize(d) + 1);
    searchRegion.SetSize(d, iLPR.GetSize(d) + jLPR.GetSize(d) - 1);
    }
  //searchRegion.SetSize(axis, 1);
  //searchRegion.SetIndex(axis, 0);
  typedef Image<bool, TImage::ImageDimension> BitmapType;
  typename BitmapType::Pointer searched = BitmapType::New();
  searched->SetRegions(searchRegion);
  searched->Allocate(true); //initialize to zero (false)

  //breadth first search starting from centroid
  std::queue<typename TImage::IndexType> uncomputed;
  uncomputed.push(ind);
  searched->SetPixel(ind, true);
  IdentifierType score, maxScore = 0;
  typename TImage::IndexType bestIndex;

  while (!uncomputed.empty())
    {
    ind = uncomputed.front();
    uncomputed.pop();
    score = Intersection(iConn, iRegionId, jConn, jRegionIds, ind);
    if (score > maxScore)
      {
      maxScore=score;
      bestIndex = ind;
      }

    //we breadth this search 
    if (!m_HeuristicAlignment || maxScore == 0 || score > maxScore / 2)
      {
      for (unsigned d = 0; d < TImage::ImageDimension; d++)
        {
        if (d == axis)
          {
          continue; //do not waste time on additional checks
          }
        ind[d] -= 1; //"left"
        if (searchRegion.IsInside(ind) && !searched->GetPixel(ind))
          {
          uncomputed.push(ind);
          searched->SetPixel(ind, true);
          }
        ind[d] += 2; //"right"
        if (searchRegion.IsInside(ind) && !searched->GetPixel(ind))
          {
          uncomputed.push(ind);
          searched->SetPixel(ind, true);
          }
        ind[d] -= 1; //return to initial
        }
      }
    }
  WriteDebug(searched, "C:\\searched.nrrd");
  return bestIndex;
}


template< typename TImage >
typename TImage::Pointer
MorphologicalContourInterpolator<TImage>
::RegionedConnectedComponents(const typename TImage::RegionType region, typename TImage::PixelType label, IdentifierType &objectCount)
{
  m_RoI->SetExtractionRegion(region);
  m_RoI->SetInput(m_Input);
  m_Binarizer->SetLowerThreshold(label);
  m_Binarizer->SetUpperThreshold(label);
  m_ConnectedComponents->Update();
  objectCount = m_ConnectedComponents->GetObjectCount();
  return m_ConnectedComponents->GetOutput();
}


template< typename TImage >
void
MorphologicalContourInterpolator<TImage>
::InterpolateBetweenTwo(int axis, TImage *out, typename TImage::IndexValueType i, typename TImage::IndexValueType j)
{
  if (i > j)
    {
    std::swap(i, j);
    }
  if (i == j || i + 1 == j)
    {
    return; //nothing to do
    }

  //compare slices i and j
  typename TImage::RegionType ri = m_TotalBoundingBox; //smaller than or equal to requested region
  ri.SetIndex(axis, i);
  ri.SetSize(axis, 1); //1 slice
  typename TImage::RegionType rj = ri;
  rj.SetIndex(axis, j);

  typename BoolImageType::Pointer eqResult = BoolImageType::New();
  typename TImage::RegionType rr = rj;
  rr.SetIndex(axis, 0);
  eqResult->CopyInformation(m_Input);
  eqResult->SetRegions(rr);
  eqResult->Allocate();

  typedef std::set<typename TImage::PixelType> LabelSetType;
  LabelSetType overlaps; //which labels have overlaps from slice i to slice j

  ImageRegionConstIterator<TImage> iti(m_Input, ri);
  ImageRegionConstIterator<TImage> itj(m_Input, rj);
  ImageRegionIterator<BoolImageType> itr(eqResult, rr);
  while (!itr.IsAtEnd())
    {
    bool eq = (iti.Get() == itj.Get() && iti.Get() != 0); //are the pixels equal and non-zero?
    itr.Set(eq);
    if (eq) //exclude background label
      {
      overlaps.insert(iti.Get());
      }
    //next pixel
    ++iti;
    ++itj;
    ++itr;
    }
  //typedef unsigned char instead of bool to write it as nrrd
  //WriteDebug<BoolImageType>(eqResult,"C:\\Temp\\eqResult.nrrd");

  //for each label with overlaps determine inter-slice region correspondences
  for (typename LabelSetType::iterator it = overlaps.begin(); it != overlaps.end(); ++it)
    {
    if (m_Label != 0 && *it != m_Label)
      continue; //label was specified, and it was not this one, so skip

    //first determine disjoint regions
    ri = m_BoundingBoxes[*it]; //even smaller than m_TotalBoundingBox
    ri.SetSize(axis, 1);
    ri.SetIndex(axis, i);
    rj = ri;
    rj.SetIndex(axis, j);
    rr = rj;
    rr.SetIndex(axis, 0);

    //execute connected components
    IdentifierType iCount, jCount;
    typename TImage::Pointer iconn = this->RegionedConnectedComponents(ri, *it, iCount);
    iconn->DisconnectPipeline();
    typename TImage::Pointer jconn = this->RegionedConnectedComponents(rj, *it, jCount);
    jconn->DisconnectPipeline();
    //WriteDebug<TImage>(iconn, "C:\\Temp\\iconn.nrrd");
    //WriteDebug<TImage>(jconn, "C:\\Temp\\jconn.nrrd");

    //go through comparison image and create correspondence pairs
    typedef std::set<std::pair<typename TImage::PixelType, typename TImage::PixelType> > PairSet;
    PairSet pairs;
    itr.SetRegion(rr);
    itr.GoToBegin();
    ImageRegionConstIterator<TImage> iti(iconn, ri);
    ImageRegionConstIterator<TImage> itj(jconn, rj);
    while (!itr.IsAtEnd())
      {
      if (itr.Get() && iti.Get() != 0 && itj.Get() != 0)
        {
        pairs.insert(std::make_pair(iti.Get(), itj.Get()));
        //std::cout << "itr:" << itr.GetIndex() <<
        //  " iti:" << iti.GetIndex() << iti.Get() <<
        //  " itj:" << itj.GetIndex() << itj.Get() << std::endl;
        }
      ++iti;
      ++itj;
      ++itr; //next pixel
      }

    typedef std::map<typename TImage::PixelType, IdentifierType> CountMap;
    CountMap iCounts, jCounts;
    typename PairSet::iterator p;
    for (p = pairs.begin(); p != pairs.end(); ++p)
      {
      iCounts[p->first]++;
      jCounts[p->second]++;
      }

    //first do extrapolation for components without overlaps
    typename CountMap::iterator iMapIt = iCounts.begin();
    for (IdentifierType ic = 1; ic <= iCount; ++ic) //component labels
      {
      if (iMapIt == iCounts.end() || ic < iMapIt->first)
        {
        Extrapolate(axis, out, *it, i, j, iconn, ic);
        }
      else //ic==iMapIt->first
        {
        ++iMapIt;
        }
      }
    typename CountMap::iterator jMapIt = jCounts.begin();
    for (IdentifierType jc = 1; jc <= jCount; ++jc) //component labels
      {
      if (jMapIt == jCounts.end() || jc < jMapIt->first)
        {
        Extrapolate(axis, out, *it, j, i, jconn, jc);
        }
      else //jc==jMapIt->first
        {
        ++jMapIt;
        }
      }

    //now handle 1 to 1 correspondences
    p = pairs.begin();
    while (p != pairs.end())
      {
      if (iCounts[p->first] == 1 && jCounts[p->second] == 1)
        {
        PixelList regionIDs;
        regionIDs.push_back(p->second);
        typename TImage::IndexType translation = Align(axis, iconn, p->first, jconn, regionIDs);
        Interpolate1to1(axis, out, *it, i, j, iconn, p->first, jconn, p->second, translation);
        iCounts.erase(p->first);
        jCounts.erase(p->second);
        pairs.erase(p++);
        }
      else
        {
        ++p;
        }
      }

    PixelList regionIDs(pairs.size()); //preallocate
    //now do 1-to-N and M-to-1 cases
    p = pairs.begin();
    while (p != pairs.end())
      {
      regionIDs.clear();

      if (iCounts[p->first] == 1) //M-to-1
        {
        for (typename PairSet::iterator rest = p; rest != pairs.end(); ++rest)
          {
          if (rest->second == p->second)
            {
            regionIDs.push_back(rest->first);
            }
          }

        typename TImage::IndexType translation = Align(axis, jconn, p->second, iconn, regionIDs);
        Interpolate1toN(axis, out, *it, j, i, jconn, p->second, iconn, regionIDs, translation);
        
        typename PairSet::iterator rest = p;
        ++rest;
        while (rest != pairs.end())
          {
          if (rest->second == p->second)
            {
            --iCounts[rest->first];
            --jCounts[rest->second];
            pairs.erase(rest++);
            }
          else
            {
            ++rest;
            }
          }
        pairs.erase(p++);
        } //M-to-1
      else if (jCounts[p->second] == 1) //1-to-N
        {
        for (typename PairSet::iterator rest = p; rest != pairs.end(); ++rest)
          {
          if (rest->first == p->first)
            {
            regionIDs.push_back(rest->second);
            }
          }

        typename TImage::IndexType translation = Align(axis, iconn, p->first, jconn, regionIDs);
        Interpolate1toN(axis, out, *it, i, j, iconn, p->first, jconn, regionIDs, translation);

        typename PairSet::iterator rest = p;
        ++rest;
        while (rest != pairs.end())
          {
          if (rest->first == p->first)
            {
            --iCounts[rest->first];
            --jCounts[rest->second];
            pairs.erase(rest++);
            }
          else
            {
            ++rest;
            }
          }
        pairs.erase(p++);
        } //1-to-N
      else
        {
        ++p;
        }
      } //1-to-N and M-to-1

    //only M-to-N correspondences remain
    //we turn each M-to-N case into m 1-to-N cases
    p = pairs.begin();
    while (p != pairs.end())
      {
      regionIDs.clear();
      for (typename PairSet::iterator rest = p; rest != pairs.end(); ++rest)
        {
        if (rest->first == p->first)
          {
          regionIDs.push_back(rest->second);
          }
        }

      typename TImage::IndexType translation = Align(axis, iconn, p->first, jconn, regionIDs);
      Interpolate1toN(axis, out, *it, i, j, iconn, p->first, jconn, regionIDs, translation);

      typename PairSet::iterator rest = p;
      ++rest;
      while (rest != pairs.end())
        {
        if (rest->first == p->first)
          {
          pairs.erase(rest++);
          }
        else
          {
          ++rest;
          }
        }
      //counts no longer matter, do not waste time deleting them
      pairs.erase(p++);
      } //M-to-N
    } //for each label with overlaps
} //void MorphologicalContourInterpolator::InterpolateBetweenTwo()


template< typename TImage >
void
MorphologicalContourInterpolator<TImage>
::InterpolateAlong(int axis, TImage *out)
{
  SliceSetType aggregate;
  if (m_Label == 0) //all labels
    {
    for (typename LabeledSlicesType::iterator it = m_LabeledSlices[axis].begin();
      it != m_LabeledSlices[axis].end(); ++it)
      {
      aggregate.insert(it->second.begin(), it->second.end());
      }
    }
  else //we only care about m_Label
    {
    aggregate = m_LabeledSlices[axis][m_Label];
    }
  typename SliceSetType::iterator prev = aggregate.begin();
  if (prev == aggregate.end())
    {
    return; //nothing to do
    }

  //set up structuring element for dilation
  typedef itk::Size<TImage::ImageDimension> SizeType;
  SizeType size;
  size.Fill(1);
  size[axis] = 0;
  m_StructuringElement.SetRadius(size);
  m_StructuringElement.CreateStructuringElement();
  m_Dilator->SetKernel(m_StructuringElement);

  typename SliceSetType::iterator it = aggregate.begin();
  for (++it; it != aggregate.end(); ++it)
    {
    InterpolateBetweenTwo(axis, out, *prev, *it);
    prev = it;
    }
}


template< typename TImage >
void
MorphologicalContourInterpolator<TImage>
::CombineInputAndInterpolate(typename TImage::Pointer interpolate)
{
  ImageRegionIterator<TImage> itO(m_Output, m_Output->GetBufferedRegion());
  ImageRegionConstIterator<TImage> itI(m_Input, m_Output->GetBufferedRegion());
  ImageRegionConstIterator<TImage> it(interpolate, m_Output->GetBufferedRegion());
  while (!it.IsAtEnd())
    {
    typename TImage::PixelType val = itI.Get();
    if (val != 0)
      {
      itO.Set(val);
      }
    else
      {
      itO.Set(it.Get());
      }

    ++it;
    ++itI;
    ++itO;
    }

  //put the output data back into the regular pipeline
  this->GraftOutput(m_Output);
  this->m_Output = ITK_NULLPTR;
}


template< typename TImage >
void
MorphologicalContourInterpolator<TImage>
::GenerateData()
{
  m_Input = TImage::New();
  m_Input->Graft(const_cast<TImage*>(this->GetInput()));
  this->AllocateOutputs();
  m_Output = TImage::New();
  m_Output->Graft(this->GetOutput());
  typename TImage::Pointer tempOut = TImage::New();
  tempOut->CopyInformation(m_Output);
  tempOut->SetRegions(m_Output->GetLargestPossibleRegion());
  
  this->DetermineSliceOrientations();

  //merge all bounding boxes
  if (m_BoundingBoxes.size() == 0)
    {
    this->GraftOutput(m_Output);
    this->m_Output = ITK_NULLPTR;
    return; //nothing to process
    }
  else
    {
    m_TotalBoundingBox = MergeBoundingBoxes(m_BoundingBoxes);
    }

  if (m_Axis == -1) //interpolate along all axes
    {
    OrientationType aggregate = OrientationType();
    aggregate.Fill(false);

    if (this->m_Label == 0)
      {
      for (typename OrientationsType::iterator it = m_Orientations.begin(); it != m_Orientations.end(); ++it)
        {
        for (unsigned int a = 0; a < TImage::ImageDimension; ++a)
          {
          aggregate[a] = aggregate[a] || it->second[a]; //any label needs interpolation along this axis
          }
        }
      }
    else
      {
      aggregate = m_Orientations[m_Label]; //we only care about this label
      }

    std::vector<typename TImage::Pointer> perAxisInterpolates;
    for (unsigned int a = 0; a < TImage::ImageDimension; ++a)
      {
      if (aggregate[a])
        {
        typename TImage::Pointer imageA = TImage::New();
        imageA->CopyInformation(m_Output);
        imageA->SetRegions(m_Output->GetRequestedRegion());
        imageA->Allocate();
        this->InterpolateAlong(a, imageA);
        perAxisInterpolates.push_back(imageA);
        }
      }

    if (perAxisInterpolates.size() == 1)
      {
      CombineInputAndInterpolate(perAxisInterpolates[0]);
      return;
      }
    //else
    std::vector<ImageRegionConstIterator<TImage> > iterators;

    for (int i = 0; i < perAxisInterpolates.size(); ++i)
      {
      ImageRegionConstIterator<TImage> it(perAxisInterpolates[i], m_Output->GetRequestedRegion());
      iterators.push_back(it);
      }

    std::vector<typename TImage::PixelType> values;
    values.reserve(perAxisInterpolates.size());

    tempOut->Allocate();
    ImageRegionIterator<TImage> it(tempOut, m_Output->GetRequestedRegion());
    while (!it.IsAtEnd())
      {
      values.clear();
      for (int i = 0; i < perAxisInterpolates.size(); ++i)
        {
        typename TImage::PixelType val = iterators[i].Get();
        if (val != 0)
          {
          values.push_back(val);
          }
        }

      if (values.size() == 0)
        {
        it.Set(0); //all were zero
        }
      else if (values.size() == 1)
        {
        it.Set(values[0]); //the only non-zero
        }
      else //median
        {
        std::nth_element(values.begin(), values.begin() + values.size() / 2, values.end());
        it.Set(values[values.size() / 2]);
        }

      //next pixel
      ++it;
      for (int i = 0; i < perAxisInterpolates.size(); ++i)
        {
        ++(iterators[i]);
        }
      }
    } //interpolate along all axes
  else //interpolate along the specified axis
    {
    tempOut->Allocate();
    this->InterpolateAlong(m_Axis, tempOut);
    }
  CombineInputAndInterpolate(tempOut);
}

}// end namespace itk

#endif //itkMorphologicalContourInterpolator_hxx
