#include "functional.h"


twoBodyFunctional::twoBodyFunctional(const potential & V,const geometrySpherical & geo) :
  _V(&V),
  _geo(geo)
{
  
};

void twoBodyFunctional::operator()(state_t & state2,const state_t & state)
  {
    stateTmp.resize(state.size());
    
    _K.evaluate(state2,state,_geo );
    _V->evaluate(stateTmp,state,_geo );
    
    state2=state2 + stateTmp;
    
  }
