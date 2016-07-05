#ifndef RLEImage_h
#define RLEImage_h

#include <utility> //std::pair
#include <vector>
#include <itkImageBase.h>
#include <itkImage.h>

/** Run-Length Encoded image.
* It saves memory for label images at the expense of processing times.
* Unsuitable for ordinary images (in which case it is counterproductive).
*
* BufferedRegion must include complete run-length lines (along X index axis).
* BufferedRegion can be smaller than LargestPossibleRegion along other axes.
*
* It is best if pixel type and counter type have the same byte size
* (for memory alignment purposes).
*
* Copied and adapted from itk::Image.
*/
template< typename TPixel, unsigned int VImageDimension = 3, typename CounterType = unsigned short >
class RLEImage : public itk::ImageBase < VImageDimension >
{

public:
    /** Standard class typedefs */
    typedef RLEImage                            Self;
    typedef itk::ImageBase < VImageDimension >  Superclass;
    typedef itk::SmartPointer< Self >           Pointer;
    typedef itk::SmartPointer< const Self >     ConstPointer;
    typedef itk::WeakPointer< const Self >      ConstWeakPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(RLEImage, ImageBase);

    /** Pixel typedef support. Used to declare pixel type in filters
    * or other operations. */
    typedef TPixel PixelType;

    /** Typedef alias for PixelType */
    typedef TPixel ValueType;

    /** First element is count of repetitions,
    * second element is the pixel value. */
    typedef std::pair<CounterType, PixelType> RLSegment;

    /** A Run-Length encoded line of pixels. */
    typedef std::vector<RLSegment> RLLine;

    /** Internal Pixel representation. Used to maintain a uniform API
    * with Image Adaptors and allow to keep a particular internal
    * representation of data while showing a different external
    * representation. */
    typedef RLLine InternalPixelType;

    //typedef PixelType IOPixelType;

    /** Dimension of the image.  This constant is used by functions that are
    * templated over image type (as opposed to being templated over pixel type
    * and dimension) when they need compile time access to the dimension of
    * the image. */
    itkStaticConstMacro(ImageDimension, unsigned int, VImageDimension);

    /** Index typedef support. An index is used to access pixel values. */
    typedef typename Superclass::IndexType      IndexType;
    typedef typename Superclass::IndexValueType IndexValueType;

    /** Offset typedef support. An offset is used to access pixel values. */
    typedef typename Superclass::OffsetType OffsetType;

    /** Size typedef support. A size is used to define region bounds. */
    typedef typename Superclass::SizeType      SizeType;
    typedef typename Superclass::SizeValueType SizeValueType;

    /** Direction typedef support. A matrix of direction cosines. */
    typedef typename Superclass::DirectionType DirectionType;

    /** Region typedef support. A region is used to specify a subset of an image.
    */
    typedef typename Superclass::RegionType RegionType;

    /** Spacing typedef support.  Spacing holds the size of a pixel.  The
    * spacing is the geometric distance between image samples. */
    typedef typename Superclass::SpacingType      SpacingType;
    typedef typename Superclass::SpacingValueType SpacingValueType;

    /** Origin typedef support.  The origin is the geometric coordinates
    * of the index (0,0). */
    typedef typename Superclass::PointType PointType;

    /** Offset typedef (relative position between indices) */
    typedef typename Superclass::OffsetValueType OffsetValueType;

    /** Allocate the image memory. The size of the image must
    * already be set, e.g. by calling SetRegions().
    * Pixel values are initialized using default constructor. */
    virtual void Allocate(bool initialize = false);

    /** Restore the data object to its initial state. This means releasing
    * memory. */
    virtual void Initialize()
    {
        // Call the superclass which should initialize the BufferedRegion ivar.
        Superclass::Initialize();
        m_OnTheFlyCleanup = true;
        myBuffer = BufferType::New();
    }

    /** Fill the image buffer with a value.  Be sure to call Allocate()
    * first. */
    void FillBuffer(const TPixel & value);

    virtual void SetLargestPossibleRegion(const RegionType & region)
    {
        Superclass::SetLargestPossibleRegion(region);
        myBuffer->SetLargestPossibleRegion(truncateRegion(region));
    }

    virtual void SetBufferedRegion(const RegionType & region)
    {
        Superclass::SetBufferedRegion(region);
        myBuffer->SetBufferedRegion(truncateRegion(region));
    }

    virtual void SetRequestedRegion(const RegionType & region)
    {
        Superclass::SetRequestedRegion(region);
        myBuffer->SetRequestedRegion(truncateRegion(region));
    }

    /** \brief Set a pixel value.
    *
    * Allocate() needs to have been called first -- for efficiency,
    * this function does not check that the image has actually been
    * allocated yet. SLOW -> Use iterators instead. */
    void SetPixel(const IndexType & index, const TPixel & value);

    /** Set a pixel value in the given line and updates segmentRemainder
    * and realIndex to refer to the same pixel.
    * Returns difference in line length which happens due to merging or splitting segments.
    * This method is used by iterators directly. */
    int SetPixel(RLLine & line, IndexValueType & segmentRemainder, IndexValueType & realIndex, const TPixel & value);

    /** \brief Get a pixel. SLOW! Better use iterators for pixel access. */
    const TPixel & GetPixel(const IndexType & index) const;

    ///** Get a reference to a pixel. Chaning it changes the whole RLE segment! */
    //TPixel & GetPixel(const IndexType & index);

    ///** \brief Access a pixel. Chaning it changes the whole RLE segment! */
    //TPixel & operator[](const IndexType & index)
    //{
    //    return this->GetPixel(index);
    //}

    /** \brief Access a pixel. This version can only be an rvalue.
    * SLOW -> Use iterators instead. */
    const TPixel & operator[](const IndexType & index) const
    {
        return this->GetPixel(index);
    }

    virtual unsigned int GetNumberOfComponentsPerPixel() const
    {
        // use the GetLength() method which works with variable length arrays,
        // to make it work with as much pixel types as possible
        PixelType p;
        return itk::NumericTraits< PixelType >::GetLength(p);
    }

    /** Typedef for the internally used buffer. */
    typedef typename itk::Image<RLLine, VImageDimension - 1> BufferType;

    /** We need to allow itk-style iterators to be constructed. */
    typename BufferType::Pointer GetBuffer() { return myBuffer; }

    /** We need to allow itk-style const iterators to be constructed. */
    typename BufferType::Pointer GetBuffer() const { return myBuffer; }

    /** Returns N-1-dimensional index, the remainder after 0-index is removed. */
    static inline typename BufferType::IndexType
        truncateIndex(const IndexType & index);

    /** Returns N-1-dimensional size, the remainder after 0-size is removed. */
    static inline typename BufferType::SizeType
        truncateSize(const SizeType & size);

    /** Returns N-1-dimensional region, the remainder after 0-index and size are removed. */
    static typename BufferType::RegionType
        truncateRegion(const RegionType & region);

    /** Merges adjacent segments with duplicate values.
    * Automatically called when turning on OnTheFlyCleanup. */
    void CleanUp() const;

    /** Should same-valued segments be merged on the fly?
    * On the fly merging usually provides better performance. */
    bool GetOnTheFlyCleanup() const { return m_OnTheFlyCleanup; }

    /** Should same-valued segments be merged on the fly?
    * On the fly merging usually provides better performance. */
    void SetOnTheFlyCleanup(bool value)
    {
        if (value == m_OnTheFlyCleanup)
            return;
        m_OnTheFlyCleanup = value;
        if (m_OnTheFlyCleanup)
            CleanUp(); //put the image into a clean state
    }


protected:
    RLEImage() : itk::ImageBase < VImageDimension >()
    {
        m_OnTheFlyCleanup = true;
        myBuffer = BufferType::New();
    }
    void PrintSelf(std::ostream & os, itk::Indent indent) const;

    virtual ~RLEImage() {}

    /** Compute helper matrices used to transform Index coordinates to
    * PhysicalPoint coordinates and back. This method is virtual and will be
    * overloaded in derived classes in order to provide backward compatibility
    * behavior in classes that did not used to take image orientation into
    * account.  */
    virtual void ComputeIndexToPhysicalPointMatrices()
    {
        this->Superclass::ComputeIndexToPhysicalPointMatrices();
    }

    /** Merges adjacent segments with duplicate values in a single line. */
    void CleanUpLine(RLLine & line) const;

private:
    bool m_OnTheFlyCleanup; //should same-valued segments be merged on the fly

    RLEImage(const Self &);          //purposely not implemented
    void operator=(const Self &); //purposely not implemented

    /** Memory for the current buffer. */
    mutable typename BufferType::Pointer myBuffer;
};


#ifndef ITK_MANUAL_INSTANTIATION
#include "RLEImage.txx"
#endif

#endif //RLEImage_h
