#ifndef utilities_h_
#define utilities_h_

// itk
#include "itkAffineTransform.h"
#include "itkImage.h"

// vnl
#include "vnl/vnl_matrix.h"


namespace kalmanAtlas
{
  /**
   * readImage
   */
  template< typename itkImage_t >
  typename itkImage_t::Pointer readImage(const char *fileName);

  /**
   * writeImage
   */
  template< typename itkImage_t > void writeImage(typename itkImage_t::Pointer img, const char *fileName);

  /**
   * Cast image pixel type
   */
  template< typename inputPixel_t, typename outputPixel_t > 
  typename itk::Image<outputPixel_t, 3 >::Pointer
  castItkImage( typename itk::Image<inputPixel_t, 3>::Pointer inputImage );


  /**
   * Read a series of images.
   */
  template< typename itkImage_t > 
  std::vector< typename itkImage_t::Pointer >
  readImageSeries( const std::vector< std::string >& imageNameList );

  /*============================================================
   * readTextLineToListOfString   
   */
  template<typename TNull>
  std::vector< std::string > readTextLineToListOfString(const char* textFileName);

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


  /**
   * Post process probabiliry map: Multiply by 200, then convert to
   * uchar image. This will make the final result to be easier
   * thresholded by Slicer.
   */
  template<typename image_t>
  typename itk::Image<unsigned char, 3>::Pointer
  postProcessProbabilityMap(typename image_t::Pointer probMap, typename image_t::PixelType c);


  /**
   * write a component of a vector image
   */
  template< typename itkVectorImage_t > 
  void 
  writeVectorImage(typename itkVectorImage_t::Pointer img, const char *fileName, int component);


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


  /**
   * read affine transformation file
   */
  template< typename CoordRepType, unsigned int NDimensions >
  typename itk::AffineTransform<CoordRepType, NDimensions>::Pointer  
  readAffineTransformFromFile(const char *fileName);


  /**
   * write affine transformation from file
   */
  template< typename CoordRepType, unsigned int NDimensions >
  void writeAffineTransformationToFile(typename itk::AffineTransform<CoordRepType, NDimensions>::Pointer trans, const char *fileName);

  /**
   * convert the affine trans to a 4x4 (3D) or 3x3 (2D) transformation
   * matrix under the homogeneous coordinate.
   */
  template< typename CoordRepType, unsigned int NDimensions >
  vnl_matrix<CoordRepType>
  affineTransformToMatrix(typename itk::AffineTransform<CoordRepType, NDimensions>::Pointer trans);


  /**
   * convert a 4x4 (3D) or 3x3 (2D) transformation matrix under the
   * homogeneous coordinate to affine trans.
   */
  template< typename CoordRepType, unsigned int NDimensions >
  typename itk::AffineTransform<CoordRepType, NDimensions>::Pointer 
  matrixToAffineTransform(vnl_matrix<CoordRepType> m, const typename itk::AffineTransform<CoordRepType, NDimensions>::ParametersType& fixedParameter);

  /**
   * Converte the homogeneous matrix to a state vector. In 3D, this
   * takes a 4x4 matrix to a 16-D vector.
   */  
  template< typename CoordRepType>
  vnl_vector<CoordRepType>
  homogeneousMatrixToStateVector(const vnl_matrix<CoordRepType>& m);


  /**
   * Converte the a state vector to a homogeneous matrix. For 3D
   * image, this takes a 16-D vector to a 4x4 matrix. 
   */  
  template< typename CoordRepType>
  vnl_matrix<CoordRepType>  
  stateVectorToHomogeneousMatrix(const vnl_vector<CoordRepType>& x);


  /**
   * Converte the homogeneous matrix to a state transition matrix. In
   * 3D, this takes a 4x4 matrix to a 16x16 matrix.
   */  
  template< typename CoordRepType>
  vnl_matrix<CoordRepType>
  homogeneousMatrixToStateTransitionMatrix(vnl_matrix<CoordRepType> m);


  /**
   * Converte a state transition matrix to a homogeneous matrix. In
   * 3D, this takes a 16x16 matrix to 4x4 matrix.
   */  
  template< typename CoordRepType>
  vnl_matrix<CoordRepType>  
  stateTransitionMatrixToHomogeneousMatrix(vnl_vector<CoordRepType> x);

  template< typename CoordRepType>
  vnl_matrix<CoordRepType>  
  invertDiagonalMatrix(const vnl_matrix<CoordRepType>& m);


  /**
   * Duplicate affine trans
   */
  template< typename CoordRepType, unsigned int NDimensions >
  typename itk::AffineTransform<CoordRepType, NDimensions>::Pointer  
  duplicateAffineTrans(typename itk::AffineTransform<CoordRepType, NDimensions>::Pointer trans);



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







}// afibReg


#include "utilities.hxx"

#endif
