#ifndef utilitiesImage_h_
#define utilitiesImage_h_

// itk
#include "itkAffineTransform.h"
#include "itkImage.h"

// vnl
#include "vnl/vnl_matrix.h"


namespace kalmanAtlas
{
  /**
   * Cast image pixel type
   */
  template< typename inputPixel_t, typename outputPixel_t > 
  typename itk::Image<outputPixel_t, 3 >::Pointer
  castItkImage( typename itk::Image<inputPixel_t, 3>::Pointer inputImage );

  template<typename image_t>
  double getVol(typename image_t::Pointer img, typename image_t::PixelType thld = 0);

  /**
   * binarilize image
   */
  template<typename input_image_t, typename output_image_t>
  typename output_image_t::Pointer
  binarilizeImage(typename input_image_t::Pointer input,                 \
                  typename input_image_t::PixelType thld,                         \
                  typename output_image_t::PixelType insideValue = 1);


  template<typename input_image_t, typename output_image_t>
  typename output_image_t::Pointer
  binarilizeImage(typename input_image_t::Pointer input,                 \
                  typename input_image_t::PixelType lowerT,                         \
                  typename input_image_t::PixelType upperT, \
                  typename output_image_t::PixelType insideValue = 1,               \
                  typename output_image_t::PixelType outsideValue = 0);


  template<typename input_image_t, typename output_image_t>
  typename output_image_t::Pointer
  thld3(typename input_image_t::Pointer input,                            \
        typename input_image_t::PixelType lowerT,                         \
        typename input_image_t::PixelType upperT, \
        typename output_image_t::PixelType insideValue = 1,               \
        typename output_image_t::PixelType outsideValue = 0);

  template<typename input_image_t, typename output_image_t>
  typename output_image_t::Pointer
  matchHistrogram(typename input_image_t::Pointer referenceImage, typename input_image_t::Pointer imageToChange);

  template<typename input_image_t, typename output_image_t>
  std::vector< typename output_image_t::Pointer >
  matchHistrogram(typename input_image_t::Pointer referenceImage,       \
                  const std::vector< typename input_image_t::Pointer >& imageListToChange);


  /**
   * Compute the non-zero region of the image
   */
  template<typename image_t>
  typename image_t::RegionType
  computeNonZeroRegion(typename image_t::Pointer img);

  /**
   * Enlarge the non-zero region so that the region is not too tightly around the non-zero reigon
   */
  template<typename image_t>
  typename image_t::RegionType
  enlargeNonZeroRegion(typename image_t::Pointer img, typename image_t::RegionType nonZeroRegion);


  /**
   * Crop the mask by its non-zero region
   */
  template<typename MaskImageType >
  typename MaskImageType::Pointer
  cropNonZeroRegionFromImage(typename MaskImageType::Pointer mask);


  /*********************************************************************************
   * About atals segmentation
   *********************************************************************************/

  /**
   * Extract the ROI from the image using the region
   */
  template<typename image_t>
  typename image_t::Pointer
  extractROI(typename image_t::Pointer img, typename image_t::RegionType region);


  /** 
   * Generate an all-zero image the same size/origin/spacing/etc. as
   * referenceImg, inside of whick, the roiRegion is the roiImg
   */
  template<typename image_t>
  typename image_t::Pointer
  antiExtractROI(typename image_t::ConstPointer roiImg,         \
                 const typename image_t::RegionType roiRegion,  \
                 typename image_t::ConstPointer referenceImg);


  template<typename image_t>
  typename image_t::RegionType
  enlargeNonZeroRegionByOnePixel(typename image_t::Pointer img, typename image_t::RegionType nonZeroRegion);


  template<typename MaskImageType>
  typename MaskImageType::Pointer
  computeSurroundingMask(typename MaskImageType::Pointer maskImage, double d);

  template<typename PriorImageType>
  typename PriorImageType::Pointer
  computeSurroundingPrior(typename PriorImageType::Pointer img);


  template<typename ImageType>
  typename ImageType::Pointer
  multiplyImageByConst(typename ImageType::Pointer img, double c);
  


}// kalmanAtlas


#include "utilitiesImage.hxx"

#endif
