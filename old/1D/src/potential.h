#include "traits.h"

class geometrySpherical;

class potential
{
public:
  virtual void evaluate( state_t & state2,const state_t & state,const geometrySpherical & geo) const =0;
private:
  
};


class squareWellPotential : potential
{
  squareWellPotential(real_t R0,real_t V0);
  
  void evaluate( state_t & state2,const state_t & state,const geometrySpherical & geo)const ;
  
private:
  
  real_t _R0,_V0;
  
  
};
