#include "kineticEnergy.h"
#include "geometry.h"


void kineticEnergy::evaluate( state_t & state2, const state_t & state,const geometrySpherical & geo)
{
  real_t invSizeStepSq=1./ (geo.sizeStep() * geo.sizeStep());
  assert(state2.size()==state.size());
  
  // compute the second derivative in the bulk
  for (int i=1;i<state2.size()-1;i++)
    {
      state2(i)= ( state(i+1) + state(i-1) - 2*state(i) )*invSizeStepSq;
    };
  // compute left second derivative
  state2(0)=( state(0) -2*state(1 ) + state(2) )*invSizeStepSq;
  
  // compute right second derivative
  state2(state2.size()-1)=( state(state.size() -1 ) -2*state(state.size()-2 ) + state( state.size()-3 ) )*invSizeStepSq;
  
}


