#include "stepper.h"


void eulerStepper::evolve( 
	MultiFab & state_new_real, MultiFab & state_new_imag,
	MultiFab & state_old_real,  MultiFab & state_old_imag,
	MLPoisson & laplacianOperatorReal, MLPoisson & laplacianOperatorImag ,
	Geometry & geom ,    const Vector<BCRec> & bc,
	Real time, Real dt)
{
	//exchange boundaries conditions
	//state_old_real.FillBoundary(geom.periodicity());
	
	//state_old_real.setDomainBndry(0, 0, 1, geom);
	//state_old_imag.setDomainBndry(0, 0, 1, geom);

    //FillDomainBoundary(state_old_real, geom, bc);
    //FillDomainBoundary(state_old_imag, geom, bc);

	evaluate(
		state_new_real,state_new_imag,
		state_old_real,state_old_imag,
		time,geom,laplacianOperatorReal,laplacianOperatorImag
		);

	state_new_real.mult(-dt);

	state_new_real.plus(state_old_real,0,1,0);
	
}
