#ifndef learnImage_hxx_
#define learnImage_hxx_

//std
#include <cstdlib>
#include <vector>

// itk
#include "itkImage.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkSubtractImageFilter.h"

// local
#include "learnImage.h"
#include "linearInterp.h"
#include "utilitiesImage.hxx"
#include "utilitiesMath.h"

namespace kalmanAtlas
{
  template<typename ImageType, typename MaskImageType>
  void computePDFUnderMask(std::vector< typename ImageType::Pointer > imgList, \
                           std::vector< typename MaskImageType::Pointer> maskList, \
                           double* pdfx, double* pdf, long pdfn, double kernelSizeFactor)
  {
    std::size_t nImg = imgList.size();
    if (nImg != maskList.size() )
      {
        std::cerr<<"Error: nImg != maskList.size()\n";
        abort();
      }

    std::vector<double> data; // store the pixel intensity under the mask

    for (std::size_t iImg = 0; iImg < nImg; ++iImg)
      {
        typename ImageType::Pointer img = imgList[iImg];
        typename MaskImageType::Pointer mask = maskList[iImg];

        typename ImageType::RegionType region = img->GetLargestPossibleRegion();

        if (region != mask->GetLargestPossibleRegion())
          {
            std::cerr<<"Error: img->GetLargestPossibleRegion() != mask->GetLargestPossibleRegion()\n";
            abort();
          }

        typedef itk::ImageRegionConstIterator<MaskImageType> MaskImageRegionConstIterator;
        MaskImageRegionConstIterator maskIt(mask, region );

        typedef itk::ImageRegionConstIterator<ImageType> ImageRegionConstIterator;
        ImageRegionConstIterator imgIt(img, region );

        maskIt.GoToBegin();
        imgIt.GoToBegin();
        for (; !imgIt.IsAtEnd(); ++imgIt, ++maskIt)
          {
            if (maskIt.Get() != 0)
              {
                data.push_back(static_cast<double>(imgIt.Get()));
              }
          }
      }

    kde_via_diffusion<double>(&data[0], data.size(),                    \
                              *min_element(data.begin(), data.end()), *max_element(data.begin(), data.end()), \
                              static_cast<std::size_t>(pdfn),           \
                              pdfx, pdf, kernelSizeFactor);

    return;
  }



  template<typename ImageType, typename MaskImageType>
  void computePDFAroundMask(std::vector< typename ImageType::Pointer > imgList, \
                            std::vector< typename MaskImageType::Pointer> maskList, \
                            double borderRadius,                        \
                            double* pdfx, double* pdf, long pdfn, double kernelSizeFactor)
  {
    std::size_t nImg = imgList.size();
    if (nImg != maskList.size() )
      {
        std::cerr<<"Error: nImg != maskList.size()\n";
        abort();
      }

    std::vector< typename MaskImageType::Pointer> surroundMaskList(nImg);
#pragma omp parallel for
    for (std::size_t iImg = 0; iImg < nImg; ++iImg)
      {
        surroundMaskList[iImg] = computeSurroundingMask<MaskImageType>(maskList[iImg], borderRadius);
      }

    computePDFUnderMask<ImageType, MaskImageType>(imgList, surroundMaskList, pdfx, pdf, pdfn, kernelSizeFactor);

    return;
  }

  template<typename ImageType, typename LikelihoodImageType>
  typename LikelihoodImageType::Pointer 
  computeLikelihoodImage(typename ImageType::Pointer img, const double* pdfx, const double* pdf, long pdfn)
  {
    typename ImageType::RegionType region = img->GetLargestPossibleRegion();

    typename LikelihoodImageType::Pointer likelihoodImg = LikelihoodImageType::New();
    likelihoodImg->SetRegions(region);
    likelihoodImg->Allocate();
    likelihoodImg->FillBuffer(0.0);
    likelihoodImg->CopyInformation(img);


    typedef itk::ImageRegionConstIterator<ImageType> ImageRegionConstIterator;
    ImageRegionConstIterator imgIt(img, region );

    typedef itk::ImageRegionIterator<LikelihoodImageType> LikelihoodImageRegionIterator;
    LikelihoodImageRegionIterator likelihoodImgIt(likelihoodImg, region);

    imgIt.GoToBegin();
    likelihoodImgIt.GoToBegin();
    for (; !imgIt.IsAtEnd(); ++imgIt, ++likelihoodImgIt)
      {
        double l = interpLinear1<double>(pdfx, pdfn, pdf, static_cast<double>(imgIt.Get()));
        likelihoodImgIt.Set(l);
      }
    
    return likelihoodImg;
  }



  /**
   * Given the pdf, compute the likelihood of a given image.
   */
  template<typename ImageType, typename LabelImageType, typename LikelihoodImageType>
  typename LikelihoodImageType::Pointer 
  computeLikelihoodImage_target(typename ImageType::Pointer img,    \
                                const std::vector< typename ImageType::Pointer >& trainingImageList, \
                                const std::vector< typename LabelImageType::Pointer >& labelImageList, \
                                long numPDFSmaple, double kernelSizeFactor)
  {
    double* pdfx_target = new double[numPDFSmaple];
    double* pdf_target = new double[numPDFSmaple];

    std::cout<<"computePDFUnderMask...\n"<<std::flush;
    computePDFUnderMask<ImageType, LabelImageType>(trainingImageList, labelImageList, pdfx_target, pdf_target, numPDFSmaple, kernelSizeFactor);

    std::cout<<"computeLikelihoodImage target...\n"<<std::flush;
    typename LikelihoodImageType::Pointer likelihoodImage_target                        \
      = computeLikelihoodImage<ImageType, LikelihoodImageType>(img, pdfx_target, pdf_target, numPDFSmaple);

    delete[] pdfx_target;
    delete[] pdf_target;

    return likelihoodImage_target;
  }

  /**
   * Given the pdf, compute the likelihood of a given image.
   */
  template<typename ImageType, typename LabelImageType, typename LikelihoodImageType>
  typename LikelihoodImageType::Pointer 
  computeLikelihoodImage_surrounding(typename ImageType::Pointer img,    \
                                     const std::vector< typename ImageType::Pointer >& trainingImageList, \
                                     const std::vector< typename LabelImageType::Pointer >& labelImageList, \
                                     double borderRadius,
                                     long numPDFSmaple, double kernelSizeFactor)
  {
    double* pdfx_surrounding = new double[numPDFSmaple];
    double* pdf_surrounding = new double[numPDFSmaple];

    //std::cout<<"computePDFAroundMask...\n"<<std::flush;
    computePDFAroundMask<ImageType, LabelImageType>(trainingImageList, labelImageList, borderRadius, \
                                                    pdfx_surrounding, pdf_surrounding, numPDFSmaple, kernelSizeFactor);

    //std::cout<<"computeLikelihoodImage surrounding...\n"<<std::flush;
    typename LikelihoodImageType::Pointer likelihoodImage_surrounding                        \
      = computeLikelihoodImage<ImageType, LikelihoodImageType>(img, pdfx_surrounding, pdf_surrounding, numPDFSmaple);

    delete[] pdfx_surrounding;
    delete[] pdf_surrounding;

    return likelihoodImage_surrounding;
  }


  /**
   * Given the pdf, compute the likelihood of a given image.
   */
  template<typename ImageType, typename LabelImageType, typename LikelihoodImageType>
  typename LikelihoodImageType::Pointer 
  computeLikelihoodDifferenceImage(typename ImageType::Pointer img,    \
                                   const std::vector< typename ImageType::Pointer >& trainingImageList, \
                                   const std::vector< typename LabelImageType::Pointer >& labelImageList, \
                                   double borderRadius,          \
                                   long numPDFSmaple, double kernelSizeFactor)
  {

    typename LikelihoodImageType::Pointer lkhdImg_target                \
      = computeLikelihoodImage_target< ImageType, LabelImageType, LikelihoodImageType>(img, \
                                                                                       trainingImageList, \
                                                                                       labelImageList, \
                                                                                       numPDFSmaple, kernelSizeFactor);

    typename LikelihoodImageType::Pointer lkhdImg_sur                \
      = computeLikelihoodImage_surrounding< ImageType, LabelImageType, LikelihoodImageType>(img, \
                                                                                            trainingImageList, \
                                                                                            labelImageList, \
                                                                                            borderRadius, \
                                                                                            numPDFSmaple, kernelSizeFactor);

    typedef itk::SubtractImageFilter<LikelihoodImageType, LikelihoodImageType > SubtractImageFilterType;
    typename SubtractImageFilterType::Pointer subtractFilter = SubtractImageFilterType::New ();
    subtractFilter->SetInput1(lkhdImg_target);
    subtractFilter->SetInput2(lkhdImg_sur);
    subtractFilter->Update();

    return subtractFilter->GetOutput();
  }



}// kalmanAtlas



#endif
