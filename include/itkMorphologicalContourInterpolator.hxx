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
#include <utility>
#include <algorithm>

namespace itk
{
template<class TImage>
typename MorphologicalContourInterpolator<TImage>::SliceSetType
MorphologicalContourInterpolator<TImage>::
  GetLabeledSliceIndices(unsigned int axis)
{
  return m_LabeledSlices[axis];
}

template<class TImage>
void MorphologicalContourInterpolator<TImage>::
  SetLabeledSliceIndices(unsigned int axis, SliceSetType indices)
{
  m_LabeledSlices[axis]=indices;
}

template<class TImage>
void MorphologicalContourInterpolator<TImage>::
SetLabeledSliceIndices(unsigned int axis, std::vector<typename TImage::IndexValueType> indices)
{
  m_LabeledSlices[axis] = SliceSetType().insert(indices.begin(), indices.end());
}

template<class TImage>
void MorphologicalContourInterpolator<TImage>::
  ExpandRegion(typename TImage::RegionType &region, typename TImage::IndexType index)
{
  for (unsigned int a = 0; a < TImage::ImageDimension; ++a)
  {
    if (region.GetIndex(a) > index[a])
    {
      region.SetIndex(a, index[a]);
    }
    else if (region.GetIndex(a) + region.GetSize(a) <= index[a])
    {
      region.SetSize(a, index[a] - region.GetIndex(a) + 1);
    }
    //else it is already within
  }
}

template<class TImage>
void MorphologicalContourInterpolator<TImage>::DetermineSliceOrientations()
{
  typename const TImage *input = this->GetInput();
  typename TImage *output = this->GetOutput();
  m_LabeledSlices.clear();
  m_LabeledSlices.resize(TImage::ImageDimension); //initialize with empty sets
  m_BoundingBoxes.clear();
  m_Orientations.clear();

  typename TImage::RegionType region = output->GetRequestedRegion();
  typename TImage::RegionType largestPossibleRegion = input->GetLargestPossibleRegion();
  for (unsigned int a = 0; a < TImage::ImageDimension; ++a)
  {
    if (region.GetIndex(a) < largestPossibleRegion.GetIndex(a) + 1)
    {
      region.SetIndex(a, largestPossibleRegion.GetIndex(a) + 1);
    }
    if (region.GetIndex(a) + region.GetSize(a) > largestPossibleRegion.GetIndex(a) + largestPossibleRegion.GetSize(a) - 1)
    {
      region.SetSize(a, largestPossibleRegion.GetIndex(a) + largestPossibleRegion.GetSize(a) - 1 - region.GetIndex(a));
    }
  }
  ImageRegionConstIteratorWithIndex<TImage> it(input, region);

  OrientationType orientations = OrientationType();
  orientations.Fill(false);

  while (!it.IsAtEnd())
  {
    typename TImage::IndexType indPrev, indNext;
    typename const TImage::IndexType ind = it.GetIndex();
    typename const TImage::PixelType val = input->GetPixel(ind);
    if (val != 0 || (m_Label != 0 && val == m_Label))
    {
      typename TImage::RegionType boundingBox1;
      boundingBox1.SetIndex(ind);
      for (unsigned int a = 0; a < TImage::ImageDimension; ++a)
        boundingBox1.SetSize(a, 1);
      std::pair<BoundingBoxesType::iterator, bool> resBB
        = m_BoundingBoxes.insert(std::make_pair(val, boundingBox1));
      if (!resBB.second) //include this index in existing BB
        ExpandRegion(resBB.first->second, ind);

      std::pair<OrientationsType::iterator, bool> res
        = m_Orientations.insert(std::make_pair(val, orientations));
      OrientationsType::iterator oRef = res.first;
      unsigned int cTrue = 0, cAdjacent = 0, axis;
      for (unsigned int a = 0; a < TImage::ImageDimension; ++a)
      {
        indPrev = ind;
        indPrev[a]--;
        indNext = ind;
        indNext[a]++;
        const typename TImage::PixelType prev = input->GetPixel(indPrev);
        const typename TImage::PixelType next = input->GetPixel(indNext);
        if (prev == 0 && next == 0)
        {
          axis = a;
          cTrue++;
        }
        else if (prev == val && next == val)
        {
          cAdjacent++;
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
    ++it;
  }

  //widen bounding boxes by 1
  region = output->GetRequestedRegion();
  for (BoundingBoxesType::iterator it = m_BoundingBoxes.begin(); it != m_BoundingBoxes.end(); ++it)
  {
    for (unsigned int a = 0; a < TImage::ImageDimension; ++a)
    {
      if (it->second.GetIndex(a) - 1 <= region.GetIndex(a))
        it->second.SetIndex(a, it->second.GetIndex(a) - 1);
      if (it->second.GetIndex(a) + it->second.GetSize(a) + 1
        <= region.GetIndex(a) + region.GetSize(a))
        it->second.SetSize(a, it->second.GetSize(a) + 1);
    }
  }
}

template< class TImage>
void MorphologicalContourInterpolator<TImage>::
  InterpolateBetweenTwo(int axis, typename TImage *out,
typename TImage::IndexValueType i, typename TImage::IndexValueType j)
{
  //determine inter-slice region correspondences
  //then do one of the three cases
  ; //for breakpoint
  //throw "todo";
}

template< class TImage>
void MorphologicalContourInterpolator<TImage>::
  InterpolateAlong(int axis, typename TImage *out)
{
  SliceSetType aggregate;
  if (m_Label == 0) //all labels
  {
    for (typename LabeledSlicesType::iterator it = m_LabeledSlices[axis].begin();
      it!= m_LabeledSlices[axis].end(); ++it)
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
  
  #pragma omp parallel
  {
    typename SliceSetType::iterator it = aggregate.begin();
    for (++it; it != aggregate.end(); ++it)
    {
      #pragma omp single nowait
      {
        InterpolateBetweenTwo(axis, out, *prev, *it);
      }
      prev = it;
    }
  }
}

template< class TImage>
void MorphologicalContourInterpolator<TImage>::GenerateData()
{
  typename TImage::ConstPointer input = this->GetInput();
  typename TImage::Pointer output = this->GetOutput();
  this->AllocateOutputs();
  output->FillBuffer(0); //clear the image now, because interpolation is optimized using bounding boxes
  
  this->DetermineSliceOrientations();

  if (m_Axis == -1)
  { 
    OrientationType aggregate = OrientationType();
    aggregate.Fill(false);

    if (this->m_Label == 0)
    {
      for (OrientationsType::iterator it = m_Orientations.begin(); it != m_Orientations.end(); ++it)
        for (unsigned int a = 0; a < TImage::ImageDimension; ++a)
          aggregate[a] = aggregate[a] || it->second[a]; //any label needs interpolation along this axis
    }
    else
      aggregate = m_Orientations[m_Label]; //we only care about this label

    std::vector<TImage::Pointer> perAxisInterpolates;
    #pragma omp parallel for
    for (unsigned int a = 0; a < TImage::ImageDimension; ++a)
    {
      if (aggregate[a])
      {
        TImage::Pointer imageA = TImage::New();
        imageA->CopyInformation(output);
        imageA->SetRegions(output->GetRequestedRegion());
        imageA->Allocate();
        this->InterpolateAlong(a, imageA);
        perAxisInterpolates.push_back(imageA);
      }
    }

    if (perAxisInterpolates.size() == 1)
    {
      output = perAxisInterpolates[0];
      return;
    }
    //else
    std::vector<ImageRegionConstIterator<TImage> > iterators;

    for (int i = 0; i < perAxisInterpolates.size(); i++)
    {
      ImageRegionConstIterator<TImage> it(perAxisInterpolates[i], output->GetRequestedRegion());
      iterators.push_back(it);
    }

    std::vector<TImage::PixelType> values;
    values.reserve(perAxisInterpolates.size());

    ImageRegionIterator<TImage> it(output, output->GetRequestedRegion());
    while (!it.IsAtEnd())
    {
      values.clear();
      for (int i = 0; i < perAxisInterpolates.size(); i++)
      {
        TImage::PixelType val = iterators[i].Value();
        if (val != 0)
          values.push_back(val);
      }

      if (values.size() == 0)
        it.Set(0); //all were zero
      else if (values.size() == 1)
        it.Set(values[0]); //the only non-zero
      else //median
      {
        std::nth_element(values.begin(), values.begin() + values.size() / 2, values.end());
        it.Set(values[values.size() / 2]);
      }

      //next pixel
      ++it;
      for (int i = 0; i < perAxisInterpolates.size(); i++)
        ++(iterators[i]);
    }
  }
  else
  {
    this->InterpolateAlong(m_Axis, output);
  }
}
}// end namespace


#endif //itkMorphologicalContourInterpolator_hxx