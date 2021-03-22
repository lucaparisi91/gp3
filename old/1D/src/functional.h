#include "potential.h"
#include "geometry.h"
#include "kineticEnergy.h"

class twoBodyFunctional
{
  twoBodyFunctional(const potential & V, const geometrySpherical & geo);
  
  void operator()(state_t & state2,const state_t & state);
private:
  
  geometrySpherical _geo;
  kineticEnergy _K;
  const potential* _V;
  state_t stateTmp;
  
  
};
