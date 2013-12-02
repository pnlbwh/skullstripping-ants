/**
 * Use the weighted multi atlas to generate a prob map

 * 20121102
 * Yi Gao
 */


//std
#include <string>
#include <limits>

// itk
#include "itkImage.h"
#include "atlasSegANTS.h"

//local
#include "utilitiesIO.h"
#include "utilitiesImage.h"

int main(int argc, char** argv)
{
  if (argc < 5)
    {
      std::cerr<<"args: rawImageName outputProbabilityImageName trainingImageListName trainingLabelImageNameList\n";
      exit(-1);
    }

  std::string rawImageName(argv[1]);
  std::string outputProbabilityImageName(argv[2]);
  std::string trainingImageListName(argv[3]); // the text file contains the absolute path of the training images
  std::string labelImageListName(argv[4]); // the text file contains the absolute path of the training label images

  const unsigned int ImageDim = 3;
  typedef float pixel_t;
  typedef itk::Image<pixel_t, ImageDim> image_t;

  typedef itk::Image<float, ImageDim> FloatImageType;

  typedef short labelPixel_t;
  typedef itk::Image<labelPixel_t, ImageDim> labelImage_t;

  /**
   * Multi-atlas registration to generate probability map.
   */
  // Read in training images' names
  std::vector< std::string > trainingImageNameList = kalmanAtlas::readTextLineToListOfString<char>(trainingImageListName.c_str());
  std::vector< image_t::Pointer > trainingImageList = kalmanAtlas::readImageSeries<image_t>(trainingImageNameList);

  // Read in raw image
  image_t::Pointer img = kalmanAtlas::readImage<image_t>(rawImageName.c_str());

  // Read in training label images' names
  std::vector< std::string > labelImageNameList = kalmanAtlas::readTextLineToListOfString<char>(labelImageListName.c_str());
  std::vector< labelImage_t::Pointer > labelImageList = kalmanAtlas::readImageSeries<labelImage_t>(labelImageNameList);

  typedef itk::Image<float, ImageDim> FloatImageType;
  //FloatImageType::Pointer probImage                                             \
    //= kalmanAtlas::atlasSegANTS_ROI_weight<image_t, labelImage_t>(img, trainingImageList, labelImageList);
  FloatImageType::Pointer probImage                                             \
    = kalmanAtlas::atlasSegANTS_ROI_weight_ants<image_t, labelImage_t>(rawImageName, trainingImageNameList, labelImageNameList);


  // multiply output image by 100, and write result
  kalmanAtlas::writeImage<FloatImageType>(kalmanAtlas::multiplyImageByConst<FloatImageType>(probImage, 100.0), outputProbabilityImageName.c_str());
  
  return 0;
}

