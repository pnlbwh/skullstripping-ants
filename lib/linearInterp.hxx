#ifndef linearInterp_hxx_
#define linearInterp_hxx_

#include <iostream>
#include <algorithm>

#include "linearInterp.h"
#include "utilitiesMath.h"

namespace kalmanAtlas
{
  template<typename FloatingPointType>
  FloatingPointType interpLinear1(const FloatingPointType* x, long nx, const FloatingPointType* y, FloatingPointType xi)
  {
    // if (!isArrayIncreasing<FloatingPointType>(x, nx))
    //   {
    //     std::cerr<<"Error: input x not increasing.\n";
    //     abort();
    //   }

    FloatingPointType yi = 0.0;

    if (xi <= x[0] || xi >= x[nx-1])
      {
        yi = 0.0;
      }
    else
      {
        long lowIdx = static_cast<long>(std::lower_bound(&x[0], &x[nx], xi) - &x[0]);
            
        if (lowIdx == 0)
          {
            yi = x[0];
          }
        else
          {
            FloatingPointType xa = x[lowIdx-1];
            FloatingPointType xb = x[lowIdx];
            FloatingPointType xd = xb - xa;

            FloatingPointType w1 = (xi - xa)/xd;
            FloatingPointType w2 = (xb - xi)/xd;
            yi = w2*y[lowIdx-1] + w1*y[lowIdx];
          }
      }

    return yi;
  }


  template<typename FloatingPointType>
  void interpLinear1(const FloatingPointType* x, long nxi, const FloatingPointType* y, \
                     FloatingPointType* xo, long nxo, FloatingPointType* yo)
  {
    if (!isArrayIncreasing<FloatingPointType>(x, nxi))
      {
        std::cerr<<"Error: input x not increasing.\n";
        abort();
      }

    for (long i = 0; i < nxo; ++i)
      {
        yo[i] = interpLinear1<FloatingPointType>(x, nxi, y, xo[i]);
      }

    return;
  }

}// kalmanAtlas


#endif
