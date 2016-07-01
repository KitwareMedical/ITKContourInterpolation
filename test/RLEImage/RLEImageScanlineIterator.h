#ifndef RLEImageScanlineIterator_h
#define RLEImageScanlineIterator_h

#include "RLEImageScanlineConstIterator.h"
#include "RLEImageIterator.h"
#include "itkImageScanlineIterator.h"

namespace itk
{
/** \class ImageScanlineIterator
* \brief A multi-dimensional iterator templated over image type that walks a
* region of pixels, scanline by scanline or in the direction of the
* fastest axis.
*/
template< typename TPixel, unsigned int VImageDimension, typename CounterType >
class ImageScanlineIterator<RLEImage<TPixel, VImageDimension, CounterType> >
    :public ImageScanlineConstIterator<RLEImage<TPixel, VImageDimension, CounterType> >
{
public:
    /** Standard class typedefs. */
    typedef ImageScanlineIterator                Self;
    typedef ImageScanlineConstIterator<RLEImage<TPixel, VImageDimension, CounterType> > Superclass;

    /** Types inherited from the Superclass */
    typedef typename Superclass::IndexType             IndexType;
    typedef typename Superclass::SizeType              SizeType;
    typedef typename Superclass::OffsetType            OffsetType;
    typedef typename Superclass::RegionType            RegionType;
    typedef typename Superclass::ImageType             ImageType;
    typedef typename Superclass::InternalPixelType     InternalPixelType;
    typedef typename Superclass::PixelType             PixelType;

    /** Default constructor. Needed since we provide a cast constructor. */
    ImageScanlineIterator()
        :ImageScanlineConstIterator< ImageType >() {}

    /** Constructor establishes an iterator to walk a particular image and a
    * particular region of that image. */
    ImageScanlineIterator(ImageType *ptr, const RegionType & region)
        :ImageScanlineConstIterator< ImageType >(ptr, region) {}

    /** Constructor that can be used to cast from an ImageIterator to an
    * ImageScanlineIterator. Many routines return an ImageIterator but for a
    * particular task, you may want an ImageScanlineIterator.  Rather than
    * provide overloaded APIs that return different types of Iterators, itk
    * returns ImageIterators and uses constructors to cast from an
    * ImageIterator to a ImageScanlineIterator. */
    ImageScanlineIterator(const ImageIterator< ImageType > & it)
        :ImageScanlineConstIterator< ImageType >(it) {}

    /** Set the pixel value */
    void Set(const PixelType & value) const
    {
        const_cast<ImageType *>(this->m_Image.GetPointer())->
            SetPixel(*const_cast<typename ImageType::RLLine *>(this->rlLine),
            this->segmentRemainder, this->realIndex, value);
    }

    ///** Return a reference to the pixel
    //* This method will provide the fastest access to pixel
    //* data, but it will NOT support ImageAdaptors. */
    //PixelType & Value(void)
    //{
    //    return myBuffer[m_Index[2]][m_Index[1]][realIndex].second;
    //}

protected:
    /** the construction from a const iterator is declared protected
    in order to enforce const correctness. */
    ImageScanlineIterator(const ImageScanlineConstIterator< ImageType > & it)
        :ImageScanlineConstIterator< ImageType >(it) {}
    Self & operator=(const ImageScanlineConstIterator< ImageType > & it)
    {
        this->ImageScanlineConstIterator< ImageType >::operator=(it);
        return *this;
    }
};
} // end namespace itk

#endif //RLEImageScanlineIterator_h