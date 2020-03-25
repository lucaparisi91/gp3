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

	if (!isImaginaryTime)
	{
		state_new_real.mult(-dt);
		state_new_imag.mult(dt);
		std::swap(state_new_real,state_new_imag);
		
	}
	else
	{
		state_new_real.mult(-dt);
		state_new_imag.mult(-dt);
	}
	state_new_real.plus(state_old_real,0,1,0);
	state_new_imag.plus(state_old_imag,0,1,0);
}

void RK4Stepper::evolve( 
	MultiFab & state_new_real, MultiFab & state_new_imag,
	MultiFab & state_old_real,  MultiFab & state_old_imag,
	MLPoisson & laplacianOperatorReal, MLPoisson & laplacianOperatorImag ,
	Geometry & geom ,    const Vector<BCRec> & bc,
	Real time, Real dt)
{
	tmp_real.define(state_new_real.boxArray(), state_new_real.DistributionMap(), nComponents, ghosts[0]);
	tmp2_real.define(state_new_real.boxArray(), state_new_real.DistributionMap(), nComponents, ghosts[0]);
	tmp_imag.define(state_new_real.boxArray(), state_new_real.DistributionMap(), nComponents, ghosts[0]);
	tmp2_imag.define(state_new_real.boxArray(), state_new_real.DistributionMap(), nComponents, ghosts[0]);


	evaluate_complex_(
		tmp_real,tmp_imag,
		state_old_real,state_old_imag,
		time,geom,laplacianOperatorReal,laplacianOperatorImag
		);

	tmp_real.mult(-dt);
	tmp_imag.mult(-dt);

	
	amrex::MultiFab::Copy(state_new_real, state_old_real, 0, 0, nComponents,ghosts);
	amrex::MultiFab::Saxpy(state_new_real,1./6.,tmp_real,0,0,nComponents,ghosts);

	amrex::MultiFab::Copy(state_new_imag, state_old_imag, 0, 0, nComponents,ghosts);
	amrex::MultiFab::Saxpy(state_new_imag,1./6.,tmp_imag,0,0,nComponents,ghosts);


	tmp_real.mult(0.5);
	tmp_real.plus(state_old_real,0,nComponents,ghosts[0]);
	tmp_imag.mult(0.5);
	tmp_imag.plus(state_old_imag,0,nComponents,ghosts[0]);

	evaluate_complex_(
		tmp2_real,tmp2_imag,
		tmp_real,tmp_imag,
		time + dt*0.5,geom,laplacianOperatorReal,laplacianOperatorImag
		);
	tmp2_real.mult(-dt);
	amrex::MultiFab::Saxpy(state_new_real,1./3.,tmp2_real,0,0,nComponents,ghosts);
	tmp2_imag.mult(-dt);
	amrex::MultiFab::Saxpy(state_new_imag,1./3.,tmp2_imag,0,0,nComponents,ghosts);

	tmp2_real.mult(0.5);
	tmp2_real.plus(state_old_real,0,nComponents,ghosts[0]);
	tmp2_imag.mult(0.5);
	tmp2_imag.plus(state_old_imag,0,nComponents,ghosts[0]);
	evaluate_complex_(
		tmp_real,tmp_imag,
		tmp2_real,tmp2_imag,
		time + dt*0.5,geom,laplacianOperatorReal,laplacianOperatorImag
		);
	tmp_real.mult(-dt);
	tmp_imag.mult(-dt);
	amrex::MultiFab::Saxpy(state_new_real,1./3.,tmp_real,0,0,nComponents,ghosts);
	amrex::MultiFab::Saxpy(state_new_imag,1./3.,tmp_imag,0,0,nComponents,ghosts);

	tmp_real.plus(state_old_real,0,nComponents,ghosts[0]);
	tmp_imag.plus(state_old_imag,0,nComponents,ghosts[0]);
	evaluate_complex_(
		tmp2_real,tmp2_imag,
		tmp_real,tmp_imag,
		time + dt,geom,laplacianOperatorReal,laplacianOperatorImag
		);
	tmp2_real.mult(-dt);
	tmp2_imag.mult(-dt);

	amrex::MultiFab::Saxpy(state_new_real,1./6.,tmp2_real,0,0,nComponents,ghosts);
	amrex::MultiFab::Saxpy(state_new_imag,1./6.,tmp2_imag,0,0,nComponents,ghosts);
}

void RK4Stepper::evaluate_complex_( 
	MultiFab & state_new_real, MultiFab & state_new_imag,
	MultiFab & state_old_real,  MultiFab & state_old_imag,
	Real time, Geometry & geom ,  MLPoisson & laplacianOperatorReal, MLPoisson & laplacianOperatorImag )
{
	evaluate(
			state_new_real,state_new_imag,
			state_old_real,state_old_imag,
			time,geom,laplacianOperatorReal,laplacianOperatorImag
			);
	if (!isImaginaryTime)
	{
		state_new_imag.mult(-1);
		std::swap(state_new_real,state_new_imag);
	}

}