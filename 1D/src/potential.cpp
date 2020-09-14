#include "potential.h"
#include "geometry.h"




squareWellPotential::squareWellPotential(real_t R0,real_t V0):
  _R0(R0),
  _V0(V0)
  {
    
  };

void squareWellPotential::evaluate( state_t & state2,const state_t & state,const geometrySpherical & geo) const
{
  assert( state.size() == state2.size() );
  for (int i=0;i<state.size();i++)
    {
      state2(i)= geo.radius(i) <= _R0 ? _V0 : 0;
    }
  
}
