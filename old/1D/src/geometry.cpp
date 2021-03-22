#include "geometry.h"


geometrySpherical::geometrySpherical(real_t rMax,size_t bins) :
  _rMax(rMax) ,
  _bins(bins)
{
  _sizeStep=_rMax/_bins;
};


real_t norm(const state_t & state, const geometrySpherical & geo)
{
  real_t norm=0;
  for (int i=0;i<state.size();i++)
    {
      norm+=4*M_PI * geo.radius(i)*geo.radius(i)* state(i)*state(i);
    }
  return norm;
};
