#ifndef utilitiesIO_h_
#define utilitiesIO_h_

// itk
#include "itkAffineTransform.h"
#include "itkImage.h"

// vnl
#include "vnl/vnl_matrix.h"


namespace kalmanAtlas
{
  /**********************************************************************************
   * readImage
   */
  template< typename itkImage_t >
  typename itkImage_t::Pointer readImage(const char *fileName);

  /************************************************************************************
   * writeImage
   */
  template< typename itkImage_t > void writeImage(typename itkImage_t::Pointer img, const char *fileName);

  /************************************************************************************
   * Read a series of images.
   */
  template< typename itkImage_t > 
  std::vector< typename itkImage_t::Pointer >
  readImageSeries( const std::vector< std::string >& imageNameList );

  /************************************************************************************
   * readTextLineToListOfString   
   */
  template<typename TNull>
  std::vector< std::string > readTextLineToListOfString(const char* textFileName);


  /************************************************************************************
   * write a component of a vector image
   */
  template< typename itkVectorImage_t > 
  void 
  writeVectorImage(typename itkVectorImage_t::Pointer img, const char *fileName, int component);


  /************************************************************************************
   * read affine transformation file
   */
  template< typename CoordRepType, unsigned int NDimensions >
  typename itk::AffineTransform<CoordRepType, NDimensions>::Pointer  
  readAffineTransformFromFile(const char *fileName);


  /************************************************************************************
   * write affine transformation from file
   */
  template< typename CoordRepType, unsigned int NDimensions >
  void writeAffineTransformationToFile(typename itk::AffineTransform<CoordRepType, NDimensions>::Pointer trans, const char *fileName);

}// kalmanAtlas


#include "utilitiesIO.hxx"

#endif
