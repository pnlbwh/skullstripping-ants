#ifndef learnImage_h_
#define learnImage_h_

//std
#include <vector>

// itk
#include "itkImage.h"


namespace kalmanAtlas
{
  /**
   * learn the pdf of the intensity under non-zero regions of the mask
   * image
   */
  template<typename ImageType, typename MaskImageType>
  void computePDFUnderMask(std::vector< typename ImageType::Pointer > imgList, \
                           std::vector< typename MaskImageType::Pointer> maskList, \
                           double* pdfx, double* pdf, long pdfn, double kernelSizeFactor);


  template<typename ImageType, typename MaskImageType>
  void computePDFAroundMask(std::vector< typename ImageType::Pointer > imgList, \
                            std::vector< typename MaskImageType::Pointer> maskList, \
                            double borderRadius,                        \
                            double* pdfx, double* pdf, long pdfn, double kernelSizeFactor);


  /**
   * Given the pdf, compute the likelihood of a given image.
   */
  template<typename ImageType, typename LikelihoodImageType>
  typename LikelihoodImageType::Pointer 
  computeLikelihoodImage(typename ImageType::Pointer img, const double* pdfx, const double* pdf, long pdfn);


  /**
   * Given the pdf, compute the likelihood of a given image.
   */
  template<typename ImageType, typename LabelImageType, typename LikelihoodImageType>
  typename LikelihoodImageType::Pointer 
  computeLikelihoodImage_target(typename ImageType::Pointer img,    \
                                const std::vector< typename ImageType::Pointer >& trainingImageList, \
                                const std::vector< typename LabelImageType::Pointer >& labelImageList, \
                                long numPDFSmaple, double kernelSizeFactor);

  /**
   * Given the pdf, compute the likelihood of a given image.
   */
  template<typename ImageType, typename LabelImageType, typename LikelihoodImageType>
  typename LikelihoodImageType::Pointer 
  computeLikelihoodImage_surrounding(typename ImageType::Pointer img,    \
                                     const std::vector< typename ImageType::Pointer >& trainingImageList, \
                                     const std::vector< typename LabelImageType::Pointer >& labelImageList, \
                                     double borderRadius,          \
                                     long numPDFSmaple, double kernelSizeFactor);

  /**
   * Given the pdf, compute the likelihood of a given image.
   */
  template<typename ImageType, typename LabelImageType, typename LikelihoodImageType>
  typename LikelihoodImageType::Pointer 
  computeLikelihoodDifferenceImage(typename ImageType::Pointer img,    \
                                   const std::vector< typename ImageType::Pointer >& trainingImageList, \
                                   const std::vector< typename LabelImageType::Pointer >& labelImageList, \
                                   double borderRadius,          \
                                   long numPDFSmaple, double kernelSizeFactor);
  
  
}// kalmanAtlas


#include "learnImage.hxx"

#endif
