#include "traits.h"
#include <cstddef>


class geometrySpherical
{
public:
  geometrySpherical(real_t rMax,size_t bins);
  
  inline real_t radius(size_t i) const { return (i+0.5)*_sizeStep;   } ;
  
  inline real_t sizeStep() const {return _sizeStep;}
  
  inline real_t radiusSphere() const {return _rMax;}
  
private:
  
  real_t _sizeStep;
  real_t _rMax;
  size_t _bins;
  
};


real_t norm(const state_t & state, const geometrySpherical & geo);
