#ifndef atlasSegANTS_h_
#define atlasSegANTS_h_

#include "itkImage.h"


namespace kalmanAtlas
{
  /*******************************************************************************
   * Use same weight for all atalses
   */
  template<typename raw_image_t, typename label_image_t>
  typename itk::Image<float, 3>::Pointer
  atlasSegANTS_ROI(typename raw_image_t::Pointer rawImg,               \
                 const typename std::vector<typename raw_image_t::Pointer>& trainingImageList, \
                 const typename std::vector<typename label_image_t::Pointer>& labelImageList);

  template<typename raw_image_t, typename label_image_t>
  typename label_image_t::Pointer
  atlasSegANTS_ROI_threshold(typename raw_image_t::Pointer rawImg,        \
                           const typename std::vector<typename raw_image_t::Pointer>& trainingImageList, \
                           const typename std::vector<typename label_image_t::Pointer>& labelImageList, \
                           double threshold);


  /*******************************************************************************
   * Give different weight for different atals
   */
  template<typename raw_image_t, typename label_image_t>
  typename itk::Image<float, 3>::Pointer
  atlasSegANTS_ROI_weight(typename raw_image_t::Pointer rawImg,               \
                        const typename std::vector<typename raw_image_t::Pointer>& trainingImageList, \
                        const typename std::vector<typename label_image_t::Pointer>& labelImageList);

  template<typename raw_image_t, typename label_image_t>
  typename itk::Image<float, 3>::Pointer
  atlasSegANTS_ROI_weight_ants(std::string rawIMageName,
                        const std::vector<std::string>& trainingImageNameList, \
                        const std::vector<std::string>& labelImageNameList);

  template<typename raw_image_t, typename label_image_t>
  typename label_image_t::Pointer
  atlasSegANTS_ROI_weight_threshold(typename raw_image_t::Pointer rawImg, \
                                  const typename std::vector<typename raw_image_t::Pointer>& trainingImageList, \
                                  const typename std::vector<typename label_image_t::Pointer>& labelImageList, \
                                  double threshold);

}


#include "atlasSegANTS.hxx"


#endif
