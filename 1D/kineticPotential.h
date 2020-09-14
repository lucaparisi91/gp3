class kineticPotentialFiniteDifference
{
  kineticPotentialFiniteDifference(real_t m);
  
  void evaluate(state_t & d2,const state_t & psi,real_t t,geometry3DSpherical & geo) // does not compute periodic boundary conditions
  {  
    for(int i=1;i<x.size()-1 ;i++)
      {
	d2(i)=-0.5*(psi(i+1) -2 *  psi(i) + psi(i-1) )/(geom.stepSize);
      }
    
  }

  evaluateBoundaries(state_t & d2,const state_t & psi,real_t t,geometry3DSpherical & geo)
  {
    
  }

private:
  const geometry_t geo;
  
}
