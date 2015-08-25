#ifndef utilitiesImage_hxx_
#define utilitiesImage_hxx_

#include <csignal>

// itk
#include "itkAffineTransform.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkCastImageFilter.h"

#include "itkDanielssonDistanceMapImageFilter.h"

#include "itkHistogramMatchingImageFilter.h"

#include "itkImage.h"
#include "itkImageDuplicator.h"

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"

//#include "itkMultiplyByConstantImageFilter.h"

#include "itkRegionOfInterestImageFilter.h"

#include "itkTransformFactoryBase.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"

#include "itkVector.h"


// vnl
#include "vnl/vnl_matrix.h"


// local
#include "utilitiesImage.h"

namespace kalmanAtlas
{
  /**
   * castItkImage
   */
  template< typename inputPixel_t, typename outputPixel_t > 
  typename itk::Image<outputPixel_t, 3 >::Pointer
  castItkImage( typename itk::Image<inputPixel_t, 3>::Pointer inputImage )
  {
    const unsigned int Dimension = 3;

    typedef itk::Image<inputPixel_t, Dimension> inputImage_t;
    typedef itk::Image<outputPixel_t, Dimension> outputImage_t;

    typedef itk::CastImageFilter< inputImage_t, outputImage_t > itkCastFilter_t;

    typename itkCastFilter_t::Pointer caster = itkCastFilter_t::New();
    caster->SetInput( inputImage );
    caster->Update();


    return caster->GetOutput();
  }


  template<typename image_t>
  double getVol(typename image_t::Pointer img, typename image_t::PixelType thld)
  {
    typedef itk::ImageRegionConstIterator<image_t> ImageRegionConstIterator_t;
    ImageRegionConstIterator_t it(img, img->GetLargestPossibleRegion() );

    double cell = (img->GetSpacing()[0])*(img->GetSpacing()[1])*(img->GetSpacing()[2]);

    double v = 0.0;

    for (it.GoToBegin(); !it.IsAtEnd(); ++it)
      {
        typename image_t::PixelType f = it.Get();

        v += (f>thld?cell:0.0);
      }

    return v;
  }



  template<typename input_image_t, typename output_image_t>
  typename output_image_t::Pointer
  thld3(typename input_image_t::Pointer input,                            \
        typename input_image_t::PixelType lowerT,                         \
        typename input_image_t::PixelType upperT, \
        typename output_image_t::PixelType insideValue,               \
        typename output_image_t::PixelType outsideValue)
  {
    /**
     * O(x) :=    I(x) \in [lowerT, upperT] ? insideValue : outsideValue
     */

    //tst
    //   std::cout<<lowerT<<std::endl;
    //   std::cout<<upperT<<std::endl;
    //tst//

    typedef itk::BinaryThresholdImageFilter<input_image_t, output_image_t> binaryThresholdImageFilter_t;

    typename binaryThresholdImageFilter_t::Pointer thlder = binaryThresholdImageFilter_t::New();
    thlder->SetInput(input);
    thlder->SetInsideValue(insideValue);
    thlder->SetOutsideValue(outsideValue);
    thlder->SetUpperThreshold(upperT);
    thlder->SetLowerThreshold(lowerT);
    thlder->Update();
  
    return thlder->GetOutput();
  }


  template<typename input_image_t, typename output_image_t>
  typename output_image_t::Pointer
  binarilizeImage(typename input_image_t::Pointer input,                 \
                  typename input_image_t::PixelType lowerT,                         \
                  typename input_image_t::PixelType upperT, \
                  typename output_image_t::PixelType insideValue,               \
                  typename output_image_t::PixelType outsideValue)
  {
    /**
     * O(x) :=    I(x) \in [lowerT, upperT] ? insideValue : outsideValue
     */

    //tst
    //   std::cout<<lowerT<<std::endl;
    //   std::cout<<upperT<<std::endl;
    //tst//

    typedef itk::BinaryThresholdImageFilter<input_image_t, output_image_t> binaryThresholdImageFilter_t;

    typename binaryThresholdImageFilter_t::Pointer thlder = binaryThresholdImageFilter_t::New();
    thlder->SetInput(input);
    thlder->SetInsideValue(insideValue);
    thlder->SetOutsideValue(outsideValue);
    thlder->SetUpperThreshold(upperT);
    thlder->SetLowerThreshold(lowerT);
    thlder->Update();
  
    return thlder->GetOutput();
  }


  template<typename input_image_t, typename output_image_t>
  typename output_image_t::Pointer
  binarilizeImage(typename input_image_t::Pointer input,                \
                  typename input_image_t::PixelType thld,               \
                  typename output_image_t::PixelType insideValue)
  {
    typename input_image_t::PixelType lowerT = thld;
    typename input_image_t::PixelType upperT = static_cast<typename input_image_t::PixelType>(1e16);
    typename output_image_t::PixelType outsideValue = 0;

    return binarilizeImage<input_image_t, output_image_t>(input, lowerT, upperT, insideValue, outsideValue);
  }


  template<typename image_t>
  typename image_t::RegionType
  computeNonZeroRegion(typename image_t::Pointer img)
  {
    /**
     * Given the img, compute the region where outside this region,
     * the image is all zero.
     *
     * The minx, y, z are initialized as sizeX, y, z; then, whenever
     * encounter an non-zero voxel, the minx, y, z are updated
     * accordingly. Similar for maxX, y, z except that they are
     * intialized to 0, 0, 0
     */
    typedef typename image_t::RegionType imageRegion_t;
    typedef typename image_t::IndexType imageIndex_t;
    typedef typename image_t::SizeType imageSize_t;


    imageRegion_t entireRegion = img->GetLargestPossibleRegion();

    long minX = entireRegion.GetSize()[0];
    long minY = entireRegion.GetSize()[1];
    long minZ = entireRegion.GetSize()[2];

    long maxX = 0;
    long maxY = 0;
    long maxZ = 0;

    //    std::cout<<"hahaha = "<<minX<<'\t'<<minY<<'\t'<<minZ<<'\t'<<maxX<<'\t'<<maxX<<'\t'<<maxX<<'\n';

    typedef itk::ImageRegionConstIteratorWithIndex< image_t > itkImageRegionConstIteratorWithIndex_t;

    itkImageRegionConstIteratorWithIndex_t it(img, entireRegion);

    char foundNonZero = 0;

    {
      imageIndex_t idx;
      for (it.GoToBegin(); !it.IsAtEnd(); ++it)
        {
          if (it.Get() != 0)
            {
              foundNonZero = 1;

              idx = it.GetIndex();

              minX = minX<idx[0]?minX:idx[0];
              minY = minY<idx[1]?minY:idx[1];
              minZ = minZ<idx[2]?minZ:idx[2];

              maxX = maxX>idx[0]?maxX:idx[0];
              maxY = maxY>idx[1]?maxY:idx[1];
              maxZ = maxZ>idx[2]?maxZ:idx[2];
            }
        }
    }

    imageRegion_t nonZeroRegion;

    if (1 == foundNonZero)
      {
        imageIndex_t startIdx;
        startIdx[0] = minX;
        startIdx[1] = minY;
        startIdx[2] = minZ;

        imageSize_t size;
        size[0] = maxX - minX;
        size[1] = maxY - minY;
        size[2] = maxZ - minZ;

        nonZeroRegion.SetSize( size );
        nonZeroRegion.SetIndex( startIdx );
      }
    else
      {
        imageIndex_t startIdx;
        startIdx[0] = 0;
        startIdx[1] = 0;
        startIdx[2] = 0;

        imageSize_t size;
        size[0] = entireRegion.GetSize()[0];
        size[1] = entireRegion.GetSize()[1];
        size[2] = entireRegion.GetSize()[2];

        nonZeroRegion.SetSize( size );
        nonZeroRegion.SetIndex( startIdx );
      }

    
    return nonZeroRegion;
  }


  /**
   * Enlarge the region by 1/5 at each end, care is taken at the
   * boundary.
   */
  template<typename image_t>
  typename image_t::RegionType
  enlargeNonZeroRegion(typename image_t::Pointer img, typename image_t::RegionType nonZeroRegion)
  {
    typedef typename image_t::RegionType imageRegion_t;
    typedef typename image_t::IndexType imageIndex_t;
    typedef typename image_t::SizeType imageSize_t;

    imageRegion_t entireRegion = img->GetLargestPossibleRegion();
    imageSize_t entireSize = entireRegion.GetSize();

    imageIndex_t start = nonZeroRegion.GetIndex();
    imageSize_t sz = nonZeroRegion.GetSize();

    start[0] = std::max(0l, static_cast<long>(start[0] - sz[0]/5));
    start[1] = std::max(0l, static_cast<long>(start[1] - sz[1]/5));
    start[2] = std::max(0l, static_cast<long>(start[2] - sz[2]/5));

    sz[0] = std::min(entireSize[0] - start[0], 7*sz[0]/5);
    sz[1] = std::min(entireSize[1] - start[1], 7*sz[1]/5);
    sz[2] = std::min(entireSize[2] - start[2], 7*sz[2]/5);

    
    /**********************************************************************************
    {
      //tst
      std::cout<<"\t\t start =    "<<start<<std::endl<<std::flush;
      std::cout<<"\t\t entireSize =    "<<entireSize<<std::endl<<std::flush;
      std::cout<<"\t\t entireSize[1] - start[1], 7*sz[1]/5   "<<entireSize[1] - start[1]<<'\t'<<7*sz[1]/5<<'\t'<<sz[1]<<std::endl<<std::flush;
      //tst//
    }
    **********************************************************************************/

    imageRegion_t largerRegion;
    largerRegion.SetSize( sz );
    largerRegion.SetIndex( start );

    return largerRegion;
  }


  /**
   * Extract the ROI from the image using the region
   */
  template<typename image_t>
  typename image_t::Pointer
  extractROI(typename image_t::Pointer img, typename image_t::RegionType region)
  {
    typedef itk::RegionOfInterestImageFilter<image_t, image_t> itkRegionOfInterestImageFilter_t;

    typename itkRegionOfInterestImageFilter_t::Pointer ROIfilter = itkRegionOfInterestImageFilter_t::New();
    ROIfilter->SetInput( img );
    ROIfilter->SetRegionOfInterest( region );
    ROIfilter->Update();

    return ROIfilter->GetOutput();
  }


  /**
   * Crop the image by the non-zero region given by the mask 
   */
  /**
   * Crop the mask by its non-zero region
   */
  template<typename MaskImageType >
  typename MaskImageType::Pointer
  cropNonZeroRegionFromImage(typename MaskImageType::Pointer mask)
  {
    typename MaskImageType::RegionType ROIRegion = computeNonZeroRegion<MaskImageType>(mask);

    typename MaskImageType::RegionType enlargedROIRegion = enlargeNonZeroRegion<MaskImageType>(mask, ROIRegion);

    typename MaskImageType::Pointer ROIMask = extractROI<MaskImageType>(mask, enlargedROIRegion);

    return ROIMask;
  }




  /** 
   * Generate an all-zero image the same size/origin/spacing/etc. as
   * referenceImg, inside of whick, the roiRegion is the roiImg
   */
  template<typename image_t>
  typename image_t::Pointer
  antiExtractROI(typename image_t::ConstPointer roiImg, const typename image_t::RegionType roiRegion, \
                 typename image_t::ConstPointer referenceImg)
  {

    /********************************************************************************
    {
      //tst
      std::cout<<"\t in antiExtractROI\n"<<std::flush;
      std::cout<<"\t roiImg.GetLargestPossibleRegion() = "<<roiImg->GetLargestPossibleRegion()<<std::endl<<std::flush;
      std::cout<<"\t roiRegion = "<<roiRegion<<std::endl<<std::flush;
      std::cout<<"\t referenceImg.GetLargestPossibleRegion() = "<<referenceImg->GetLargestPossibleRegion()<<std::endl<<std::flush;
      //tst//
    }
    ********************************************************************************/

    typedef typename image_t::Pointer imagePointer_t;

    imagePointer_t largeImage = image_t::New();
    largeImage->SetRegions( referenceImg->GetLargestPossibleRegion() );
    largeImage->Allocate();

    largeImage->FillBuffer(0);
    largeImage->CopyInformation(referenceImg);


    typedef itk::ImageRegionIterator< image_t > itkImageRegionIterator_t;
    typedef itk::ImageRegionConstIterator< image_t > itkImageRegionConstIterator_t;

    {
      itkImageRegionConstIterator_t itROI(roiImg, roiImg->GetLargestPossibleRegion());
      itkImageRegionIterator_t itNew(largeImage, roiRegion);

      itROI.GoToBegin();
      itNew.GoToBegin();
      for (; !itROI.IsAtEnd(); ++itROI, ++itNew)
        {
          itNew.Set(itROI.Get());
        }
    }

    return largeImage;
  }



  template<typename image_t>
  typename image_t::RegionType
  enlargeNonZeroRegionByOnePixel(typename image_t::Pointer img, typename image_t::RegionType nonZeroRegion)
  {
    typedef typename image_t::RegionType imageRegion_t;
    typedef typename image_t::IndexType imageIndex_t;
    typedef typename image_t::SizeType imageSize_t;

    imageRegion_t entireRegion = img->GetLargestPossibleRegion();
    imageSize_t entireSize = entireRegion.GetSize();

    imageIndex_t start = nonZeroRegion.GetIndex();
    imageSize_t sz = nonZeroRegion.GetSize();


    start[0] = std::max(0l, static_cast<long>(start[0] - 1));
    start[1] = std::max(0l, static_cast<long>(start[1] - 1));
    start[2] = std::max(0l, static_cast<long>(start[2] - 1));

    sz[0] = std::min(entireSize[0] - start[0], sz[0] + 2);
    sz[1] = std::min(entireSize[1] - start[1], sz[1] + 2);
    sz[2] = std::min(entireSize[2] - start[2], sz[2] + 2);


    /**********************************************************************************    
    {
      //tst
      std::cout<<"\t\t start =    "<<start<<std::endl<<std::flush;
      std::cout<<"\t\t entireSize =    "<<entireSize<<std::endl<<std::flush;
      std::cout<<"\t\t entireSize[1] - start[1], 7*sz[1]/5   "<<entireSize[1] - start[1]<<'\t'<<7*sz[1]/5<<'\t'<<sz[1]<<std::endl<<std::flush;
      //tst//
    }
    ********************************************************************************/

    imageRegion_t largerRegion;
    largerRegion.SetSize( sz );
    largerRegion.SetIndex( start );

    return largerRegion;
  }


  template<typename input_image_t, typename output_image_t>
  typename output_image_t::Pointer
  matchHistrogram(typename input_image_t::Pointer referenceImage, typename input_image_t::Pointer imageToChange)
  {
    // Histogram match the images
    typedef itk::HistogramMatchingImageFilter<input_image_t, output_image_t> HMImageFilterType;
    typename HMImageFilterType::Pointer HMImageFilter = HMImageFilterType::New();

    HMImageFilter->SetReferenceImage( referenceImage );
    HMImageFilter->SetInput( imageToChange );
    HMImageFilter->SetNumberOfHistogramLevels( 100 );
    HMImageFilter->SetNumberOfMatchPoints( 15);
    HMImageFilter->ThresholdAtMeanIntensityOn();
    HMImageFilter->Update();

    return HMImageFilter->GetOutput();    
  }


  template<typename input_image_t, typename output_image_t>
  std::vector< typename output_image_t::Pointer >
  matchHistrogram(typename input_image_t::Pointer referenceImage,       \
                  const std::vector< typename input_image_t::Pointer >& imageListToChange)
  {
    long n = static_cast<long>(imageListToChange.size());

    std::vector< typename output_image_t::Pointer > outputImageList(n);

#pragma omp parallel for
    for (long i = 0; i < n; ++i)
      {
        outputImageList[i] = matchHistrogram<input_image_t, output_image_t>(referenceImage, imageListToChange[i]);
      }

    return outputImageList;
  }

  template<typename MaskImageType>
  typename MaskImageType::Pointer
  computeSurroundingMask(typename MaskImageType::Pointer maskImage, double d)
  {
    if (d <= 0 )
      {
        std::cerr<<"Error: d <= 0\n";
        abort();
      }

    /** Extract the non-zero region to form a 0-1 mask image */
    typedef itk::ImageDuplicator< MaskImageType > DuplicatorType;
    typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
    duplicator->SetInputImage(maskImage);
    duplicator->Update();

    typedef itk::ImageRegionIterator<MaskImageType> MaskImageIterator;
    typename MaskImageType::Pointer binaryImage = duplicator->GetOutput();
    MaskImageIterator itMask(binaryImage, binaryImage->GetLargestPossibleRegion());

    for (itMask.GoToBegin(); !itMask.IsAtEnd(); ++itMask)
      {
        if (itMask.Get() != 0)
          {
            itMask.Set(1);
          }
      }

    /* build distrance map from the non-zero region of the mask image */
    typedef itk::Image<float, MaskImageType::ImageDimension> FloatImageType;

    typedef itk::DanielssonDistanceMapImageFilter< MaskImageType, FloatImageType >  FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    filter->SetInput( binaryImage );
    filter->InputIsBinaryOn();
    filter->UseImageSpacingOn();
    filter->Update();

    typename FloatImageType::Pointer distMap = filter->GetOutput();

    /* Extract the bouder region */
    typedef itk::ImageRegionIterator<FloatImageType> FloatImageIterator;
    FloatImageIterator itDist(distMap, distMap->GetLargestPossibleRegion());


    itDist.GoToBegin();
    itMask.GoToBegin();
    for (; !itMask.IsAtEnd(); ++itDist, ++itMask)
      {
        if (itMask.Get() != 1 && itDist.Get() <= d)
          {
            itMask.Set(1);
          }
        else
          {
            itMask.Set(0);
          }
      }

    return binaryImage;
  }


  template<typename PriorImageType>
  typename PriorImageType::Pointer
  computeSurroundingPrior(typename PriorImageType::Pointer img)
  {
    typename PriorImageType::RegionType region = img->GetLargestPossibleRegion();
    typename PriorImageType::Pointer imgNew = PriorImageType::New();
    imgNew->SetRegions(region);
    imgNew->Allocate();
    imgNew->FillBuffer(0.0);
    imgNew->CopyInformation(img);

    
    typedef itk::ImageRegionIterator<PriorImageType> ImageRegionIterator;
    ImageRegionIterator itImg(img, region);
    ImageRegionIterator itImgNew(imgNew, region);

    itImg.GoToBegin();
    itImgNew.GoToBegin();

    for (; !itImg.IsAtEnd(); ++itImg, ++itImgNew)
      {
        itImgNew.Set(1.0 - itImg.Get());
      }
    
    return imgNew;
  }

  template<typename ImageType>
  typename ImageType::Pointer
  multiplyImageByConst(typename ImageType::Pointer img, double c)
  {
    typedef typename ImageType::Pointer ImagePointer;
    ImagePointer newImg = ImageType::New();
    newImg->SetRegions(img->GetLargestPossibleRegion());
    newImg->Allocate();
    newImg->CopyInformation(img);

    typedef itk::ImageRegionIterator< ImageType > itkImageRegionIterator_t;
    typedef itk::ImageRegionConstIterator< ImageType > itkImageRegionConstIterator_t;

    {
      itkImageRegionConstIterator_t it(img, img->GetLargestPossibleRegion());
      itkImageRegionIterator_t itNew(newImg, newImg->GetLargestPossibleRegion());

      it.GoToBegin();
      itNew.GoToBegin();
      for (; !it.IsAtEnd(); ++it, ++itNew)
        {
          itNew.Set(c*(it.Get()));
        }
    }
    
    return newImg;
  }


}// kalmanAtlas

#endif
