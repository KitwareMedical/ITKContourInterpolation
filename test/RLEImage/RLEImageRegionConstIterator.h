#ifndef RLEImageRegionConstIterator_h
#define RLEImageRegionConstIterator_h

#include "RLEImageConstIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionConstIteratorWithOnlyIndex.h"

class MultiLabelMeshPipeline;

namespace itk
{
/** \class ImageRegionConstIterator
 * \brief A multi-dimensional iterator templated over image type that walks a
 * region of pixels.
 *
 * ImageRegionConstIterator provides read-only access to image data.  It is the
 * base class for the read/write access ImageRegionIterator.
 *
 */
template< typename TPixel, unsigned int VImageDimension, typename CounterType >
class ImageRegionConstIterator<RLEImage<TPixel, VImageDimension, CounterType> >
    :public ImageConstIterator<RLEImage<TPixel, VImageDimension, CounterType> >
{
    friend class ::MultiLabelMeshPipeline;
public:
  /** Standard class typedef. */
  typedef ImageRegionConstIterator<RLEImage<TPixel, VImageDimension, CounterType> >     Self;
  typedef ImageConstIterator<RLEImage<TPixel, VImageDimension, CounterType> > Superclass;

  /** Dimension of the image that the iterator walks.  This constant is needed so
   * functions that are templated over image iterator type (as opposed to
   * being templated over pixel type and dimension) can have compile time
   * access to the dimension of the image that the iterator walks. */
  itkStaticConstMacro(ImageIteratorDimension, unsigned int, VImageDimension);

  /**
   * Index typedef support. While these were already typdef'ed in the superclass,
   * they need to be redone here for this subclass to compile properly with gcc.
   */
  /** Types inherited from the Superclass */
  typedef typename Superclass::IndexType             IndexType;
  typedef typename Superclass::SizeType              SizeType;
  typedef typename Superclass::OffsetType            OffsetType;
  typedef typename Superclass::RegionType            RegionType;
  typedef typename Superclass::ImageType             ImageType;
  typedef typename Superclass::InternalPixelType     InternalPixelType;
  typedef typename Superclass::PixelType             PixelType;

  /** Run-time type information (and related methods). */
  itkTypeMacro(ImageRegionConstIterator, ImageConstIterator);

  /** Default constructor. Needed since we provide a cast constructor. */
  ImageRegionConstIterator() :ImageConstIterator< ImageType >(){ }

  /** Constructor establishes an iterator to walk a particular image and a
   * particular region of that image. */
  ImageRegionConstIterator(const ImageType *ptr, const RegionType & region):
    ImageConstIterator< ImageType >(ptr, region) { }

  /** Constructor that can be used to cast from an ImageIterator to an
   * ImageRegionConstIterator. Many routines return an ImageIterator, but for a
   * particular task, you may want an ImageRegionConstIterator.  Rather than
   * provide overloaded APIs that return different types of Iterators, itk
   * returns ImageIterators and uses constructors to cast from an
   * ImageIterator to a ImageRegionConstIterator. */
  ImageRegionConstIterator(const ImageIterator< ImageType > & it)
  {
    this->ImageConstIterator< ImageType >::operator=(it);
    //this->ImageConstIterator< ImageType >::operator=(static_cast<const ImageConstIterator<ImageType> >(it));
  }

  /** Constructor that can be used to cast from an ImageConstIterator to an
   * ImageRegionConstIterator. Many routines return an ImageIterator, but for a
   * particular task, you may want an ImageRegionConstIterator.  Rather than
   * provide overloaded APIs that return different types of Iterators, itk
   * returns ImageIterators and uses constructors to cast from an
   * ImageIterator to a ImageRegionConstIterator. */
  ImageRegionConstIterator(const ImageConstIterator< ImageType > & it)
  {
    this->ImageConstIterator< ImageType >::operator=(it);
  }

  /** Increment (prefix) the fastest moving dimension of the iterator's index.
   * This operator will constrain the iterator within the region (i.e. the
   * iterator will automatically wrap from the end of the row of the region
   * to the beginning of the next row of the region) up until the iterator
   * tries to moves past the last pixel of the region.  Here, the iterator
   * will be set to be one pixel past the end of the region.
   * \sa operator++(int) */
  Self &
  operator++()
  {
    this->m_Index0++;

    if (this->m_Index0 >= this->m_EndIndex0)
    {
        ++(this->bi);
        if (!this->bi.IsAtEnd())
            this->SetIndexInternal(this->m_BeginIndex0);
        else
            this->m_Index0 = this->m_BeginIndex0;
        return *this;
    }

    this->segmentRemainder--;
    if (this->segmentRemainder > 0)
        return *this;

    this->realIndex++;
    this->segmentRemainder = (*this->rlLine)[this->realIndex].first;
    return *this;
  }

  /** Decrement (prefix) the fastest moving dimension of the iterator's index.
   * This operator will constrain the iterator within the region (i.e. the
   * iterator will automatically wrap from the beginning of the row of the region
   * to the end of the next row of the region) up until the iterator
   * tries to moves past the first pixel of the region.  Here, the iterator
   * will be set to be one pixel past the beginning of the region.
   * \sa operator--(int) */
  Self & operator--()
  {
      this->m_Index0--;

      if (this->m_Index0 < this->m_BeginIndex0)
      {
          --(this->bi);
          this->SetIndexInternal(this->m_EndIndex0 - 1);
          return *this;
      }

      this->segmentRemainder++;
      if (this->segmentRemainder <= (*this->rlLine)[this->realIndex].first)
          return *this;

      this->realIndex--;
      this->segmentRemainder = 1;
      return *this;
  }
};

template< typename TPixel, unsigned int VImageDimension, typename CounterType >
class ImageRegionConstIteratorWithIndex<RLEImage<TPixel, VImageDimension, CounterType> >
    :public ImageRegionConstIterator < RLEImage<TPixel, VImageDimension, CounterType> >
{
public:
    typedef RLEImage<TPixel, VImageDimension, CounterType> ImageType;

    typedef typename itk::ImageConstIterator<RLEImage<TPixel, VImageDimension, CounterType> >::RegionType RegionType;

    /** Default constructor. Needed since we provide a cast constructor. */
    ImageRegionConstIteratorWithIndex() :ImageRegionConstIterator< ImageType >(){ }

    /** Constructor establishes an iterator to walk a particular image and a
    * particular region of that image. */
    ImageRegionConstIteratorWithIndex(const ImageType *ptr, const RegionType & region) :
        ImageRegionConstIterator< ImageType >(ptr, region) { }

    void GoToReverseBegin()
    {
        this->bi.GoToEnd(); //after last pixel
        --(this->bi); //go to last valid pixel
        this->m_Index0 = this->m_EndIndex0 - 1;
        this->SetIndexInternal(this->m_Index0); //valid index required
    }

    bool IsAtReverseEnd()
    {
        return (this->m_Index0 == this->m_BeginIndex0) && this->bi.IsAtBegin();
    }

    /** Constructor that can be used to cast from an ImageIterator to an
    * ImageRegionConstIterator. Many routines return an ImageIterator, but for a
    * particular task, you may want an ImageRegionConstIterator.  Rather than
    * provide overloaded APIs that return different types of Iterators, itk
    * returns ImageIterators and uses constructors to cast from an
    * ImageIterator to a ImageRegionConstIterator. */
    ImageRegionConstIteratorWithIndex(const ImageIterator< ImageType > & it)
    {
        this->ImageRegionConstIterator< ImageType >::operator=(it);
    }

}; //no additional implementation required

template< typename TPixel, unsigned int VImageDimension, typename CounterType >
class ImageRegionConstIteratorWithOnlyIndex<RLEImage<TPixel, VImageDimension, CounterType> >
    :public ImageRegionConstIteratorWithIndex < RLEImage<TPixel, VImageDimension, CounterType> >
{
    //just inherit constructors
public:

    typedef RLEImage<TPixel, VImageDimension, CounterType> ImageType;

    typedef typename itk::ImageConstIterator<RLEImage<TPixel, VImageDimension, CounterType> >::RegionType RegionType;

    /** Default constructor. Needed since we provide a cast constructor. */
    ImageRegionConstIteratorWithOnlyIndex() :ImageRegionConstIterator< ImageType >(){ }

    /** Constructor establishes an iterator to walk a particular image and a
    * particular region of that image. */
    ImageRegionConstIteratorWithOnlyIndex(const ImageType *ptr, const RegionType & region) :
        ImageRegionConstIteratorWithIndex< ImageType >(ptr, region) { }

    /** Constructor that can be used to cast from an ImageIterator to an
    * ImageRegionConstIterator. Many routines return an ImageIterator, but for a
    * particular task, you may want an ImageRegionConstIterator.  Rather than
    * provide overloaded APIs that return different types of Iterators, itk
    * returns ImageIterators and uses constructors to cast from an
    * ImageIterator to a ImageRegionConstIterator. */
    ImageRegionConstIteratorWithOnlyIndex(const ImageIterator< ImageType > & it)
    {
        this->ImageRegionConstIterator< ImageType >::operator=(it);
    }

}; //no additional implementation required

} // end namespace itk

#endif //RLEImageRegionConstIterator_h
