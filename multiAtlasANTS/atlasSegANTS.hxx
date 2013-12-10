#ifndef atlasSegANTS_txx_
#define atlasSegANTS_txx_


#include <cstdlib>
#include <cmath>
#include <algorithm> // for max and min and copy
#include <iterator> // for ostream_iterator
#include <vector>

// itk
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkCastImageFilter.h"

// ANTS
#include "ants.h"
#include "itkAvantsMutualInformationRegistrationFunction.h"

// local
#include "atlasSegANTS.h"

#include "utilitiesImage.h"
#include "utilitiesIO.h"


#include "cmakeConfig.h"
#include <libgen.h>  // for basename

namespace kalmanAtlas
{


  template<typename raw_image_t, typename label_image_t>
  typename itk::Image<float, 3>::Pointer
  atlasSegANTS_ROI_weight_ants(std::string rawImageName,
                        const std::vector<std::string>& trainingImageNameList, \
                        const std::vector<std::string>& labelImageNameList)
  {

      const unsigned int ImageDim = 3;
      typedef float pixel_t;

      typedef float internal_pixel_t;
      typedef itk::Image<internal_pixel_t, 3> internal_image_t;

      typedef itk::Image<pixel_t, ImageDim> image_t;
      image_t::Pointer rawImg = kalmanAtlas::readImage<image_t>(rawImageName.c_str());

      unsigned long nTrainingImg = trainingImageNameList.size();
      std::vector< typename internal_image_t::Pointer > transformedLabelImages(nTrainingImg);
      std::vector<double> registrationFinalCostFunctionalValues(nTrainingImg);

      // Strip directory path (if any) and extension
      std::string rawImageNameBase = std::string(basename(&rawImageName[0]));
      int lastindex = rawImageNameBase.find_last_of("."); 
      rawImageNameBase = rawImageNameBase.substr(0, lastindex); 
      std::string mkdtemp_template = "/tmp/" + rawImageNameBase + "_XXXXXX"; 
      char* tmpdir = mkdtemp(&mkdtemp_template[0]);

      for (long iTrainingImg = 0; iTrainingImg < (long)nTrainingImg; ++iTrainingImg)
      {
          std::cout << "Image " << iTrainingImg << std::endl << std::flush;
          std::cout << "=======" << std::endl << std::endl << std::flush;
          // Non-linear registration, including affine
          std::ostringstream prefix;
          prefix << tmpdir << "/" << rawImageNameBase << "_ants" << iTrainingImg;
          std::cout << "--- ANTS: Running " << prefix.str() << "Affine.txt," << prefix.str() << "Warp.nii.gz <- " << std::endl;
          std::vector<std::string> args;
          args.push_back("3");
          args.push_back("-m");
          std::string cost = "CC[";
          cost += rawImageName;
          cost += ",";
          cost += trainingImageNameList[iTrainingImg]; //moving
          cost += ",1,5]";
          args.push_back(cost);
          args.push_back("-i");
          args.push_back("50x20x10"); //args.push_back("50x20x10");
          args.push_back("-r");
          args.push_back("Gauss[3,0]");
          args.push_back("-o");
          args.push_back(prefix.str());
          for( int i = 0; i < args.size(); ++i)
              std::cout << args[i] << " " <<std::flush;
          std::cout << std::endl;
          ants::ANTS(args, &std::cout);

          // Transform images
          std::vector<std::string> args2;
          args2.push_back("3");
          args2.push_back(trainingImageNameList[iTrainingImg]); //moving
          args2.push_back(prefix.str() + ".nii.gz");  //out image
          args2.push_back("-R");
          args2.push_back(rawImageName);
          args2.push_back(prefix.str() + "Warp.nii.gz");
          args2.push_back(prefix.str() + "Affine.txt");
          std::cout << std::endl << "--- WarpImageMultiTransform: Running " << prefix.str() << ".nii.gz <- " << std::endl;
          for( int i = 0; i < args2.size(); ++i)
              std::cout << args2[i] << " " << std::flush;
          std::cout << std::endl;
          ants::WarpImageMultiTransform(args2, &std::cout);
          args2[1] = labelImageNameList[iTrainingImg];
          args2[2] = prefix.str() + "label" + ".nii.gz";
          std::cout << "--- WarpImageMultiTransform: Running " << args2[2] << " <- " << std::endl;
          for( int i = 0; i < args2.size(); ++i)
              std::cout << args2[i] << " " << std::flush;
          ants::WarpImageMultiTransform(args2, &std::cout);
          transformedLabelImages[iTrainingImg] = kalmanAtlas::readImage<internal_image_t>(args2[2].c_str());
          std::cout << std::endl << "--- WarpImageMultiTransform: Done"  << std::endl;
          std::cout << std::endl;

          // Measure image similarity
          std::cout << "--- MeasureImageSimilarity: Running (ANTS similarity) <- " << prefix.str() << ".nii.gz," << rawImageName << std::endl << std::flush;

          //std::vector<std::string> args3;
          //args3.push_back("3");
          //args3.push_back("2");  // 2 for mutual information
          //args3.push_back(trainingImageNameList[iTrainingImg]); 
          //args3.push_back(rawImageName);
          //std::stringstream ants_output;
          //std::cout << "MeasureImageSimilarity: calling with args: "  << std::endl;
          //copy(args3.begin(), args3.end(), std::ostream_iterator<std::string>(std::cout, " "));
          //std::cout << std::endl << std::flush;
          //ants::MeasureImageSimilarity(args3, &ants_output);
          //std::cout << "MeasureImageSimilarity: done "  << std::endl << std::flush;
          //std::vector<std::string> output;
          //std::string token;
          //std::cout << "tokens: ";
          //while (ants_output >> token) 
          //{
              //output.push_back(token);
              //std::cout << token << "|" << std::flush;
          //}
          //std::cout << std::endl;
          //registrationFinalCostFunctionalValues[iTrainingImg] = ::atof(output[5].c_str());
          
          image_t::Pointer trainingImg = kalmanAtlas::readImage<image_t>(trainingImageNameList[iTrainingImg].c_str()); 

          typedef internal_image_t FixedImageType;
          typedef internal_image_t MovingImageType;
          typedef itk::Vector<float, 3>                     VectorType;
          typedef itk::Image<VectorType, 3>                 FieldType;
          typedef FieldType DisplacementFieldType;
          typedef itk::AvantsMutualInformationRegistrationFunction<FixedImageType, MovingImageType,
                                                           DisplacementFieldType>  MIMetricType;

          typename MIMetricType::RadiusType miradius;
          miradius.Fill(0);
          typename MIMetricType::Pointer mimet = MIMetricType::New();
          mimet->SetFixedImage(rawImg);
          mimet->SetMovingImage(trainingImg);
          mimet->SetRadius(miradius);
          mimet->SetGradientStep(1.e2);
          mimet->SetNormalizeGradient(false);

          double      metricvalue = 0;
          mimet->InitializeIteration();
          metricvalue = mimet->ComputeMutualInformation();

         registrationFinalCostFunctionalValues[iTrainingImg] = metricvalue;
         std::cout << "MeasureImageSimilarity: cost function: " << registrationFinalCostFunctionalValues[iTrainingImg] << std::endl;
         std::cout << "--- MeasureImageSimilarity: Done" << std::endl << std::flush;
         std::cout << std::endl;
          
      }

      std::cout<<"Finished ANTS\n"<<std::flush;


      // 2.
      typename internal_image_t::Pointer resultLabelImg = internal_image_t::New();
      resultLabelImg->SetRegions(rawImg->GetLargestPossibleRegion());
      resultLabelImg->Allocate();
      resultLabelImg->SetSpacing(rawImg->GetSpacing() );
      resultLabelImg->CopyInformation(rawImg );
      resultLabelImg->FillBuffer(0);

      typedef itk::ImageRegionConstIterator<internal_image_t> ImageRegionConstIterator_t;
      typedef itk::ImageRegionIterator<internal_image_t> ImageRegionIterator_t;



      /**
       * Fusion
       *
       * 1. from final cost functional values, compute the weights
       * 2. weighted sum of transformed label images
       */
      {
          // 1.
          double maxCostFnalValue = *max_element(registrationFinalCostFunctionalValues.begin(), \
                  registrationFinalCostFunctionalValues.end()); // worst reg

          double minCostFnalValue = *min_element(registrationFinalCostFunctionalValues.begin(), \
                  registrationFinalCostFunctionalValues.end());  // best reg

          double costFnalRange = maxCostFnalValue - minCostFnalValue;

          //     // compute factor s.t. exp(-factor*costFnalRange) = 0.1
          //     double minus_ln01 = 2.302585093; // -ln(0.1)

          // compute factor s.t. exp(-factor*costFnalRange) = 0.2

          /**
           * for the HUVA data set, 0.2 is better than 0.1, 0.01, 0.001, and not-weight (same weight average)
           */
          double minus_ln01 = -log(0.2);

          double factor = minus_ln01/costFnalRange;


          assert(nTrainingImg == registrationFinalCostFunctionalValues.size());

          std::vector<double> weights(nTrainingImg);

          for (int i = 0; i < static_cast<int>(nTrainingImg); ++i)
          {
              weights[i] = exp(factor*(minCostFnalValue - registrationFinalCostFunctionalValues[i]));
          }

          // normalize weights s.t. sum to 1
          double s = 0.0;
          for (int i = 0; i < static_cast<int>(nTrainingImg); ++i)
          {
              s += weights[i];
          }

          for (int i = 0; i < static_cast<int>(nTrainingImg); ++i)
          {
              weights[i] /= s;
              std::cout<<weights[i]<<"  ";
          }
          std::cout<<std::endl;


          // 2.
          {
              ImageRegionIterator_t resultLabelImgIter(resultLabelImg, resultLabelImg->GetLargestPossibleRegion() );

              for (int iTrainingImg = 0; iTrainingImg < static_cast<int>(nTrainingImg); ++iTrainingImg)
              {
                  resultLabelImgIter.GoToBegin();

                  typename internal_image_t::Pointer regLabelImg = transformedLabelImages[iTrainingImg];
                  ImageRegionConstIterator_t regLabelImgIter(regLabelImg, regLabelImg->GetLargestPossibleRegion() );
                  regLabelImgIter.GoToBegin();

                  for (; !regLabelImgIter.IsAtEnd(); ++regLabelImgIter, ++resultLabelImgIter)
                  {
                      resultLabelImgIter.Set(resultLabelImgIter.Get() + weights[iTrainingImg]*regLabelImgIter.Get() );
                  }
              }
          } // 2.//

      } // Fusion


    std::cout << "Finished atlasSegANTS_ROI_weight_ants" << std::endl;

      // cast to float and output
      return kalmanAtlas::castItkImage<internal_pixel_t, float>( resultLabelImg ) ;
  }


}// namespace
#endif
