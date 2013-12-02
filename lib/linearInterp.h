#ifndef linearInterp_h_
#define linearInterp_h_

namespace kalmanAtlas
{
  /**
   * Linear interpolation in 1D. Note: for efficiency, this does NOT
   * check strict-increasing of the x
   */
  template<typename FloatingPointType>
  FloatingPointType interpLinear1(const FloatingPointType* x, long nx, const FloatingPointType* y, FloatingPointType xi);

  /**
   * Linear interpolation in 1D. Note: this checks strict-increasing
   * of the x once.
   */
  template<typename FloatingPointType>
  void interpLinear1(const FloatingPointType* x, long nxi, const FloatingPointType* y, \
                     FloatingPointType* xo, long nxo, FloatingPointType* yo);

}// kalmanAtlas


#include "linearInterp.hxx"

#endif
