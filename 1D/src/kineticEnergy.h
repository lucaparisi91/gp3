#ifndef KINETIC_ENERGY_H
#define KINETIC_ENERGY_H

#include "traits.h"

class geometrySpherical;





class kineticEnergy
{
public:
  kineticEnergy(){};
  
  void evaluate( state_t & state2, const state_t & state,const geometrySpherical & geo);
  
};


#endif
