#include "stepper.h"

void eulerStepper::evolve( 
	MultiFab & state_new_real, MultiFab & state_new_imag,
	MultiFab & state_old_real,  MultiFab & state_old_imag,
	Real time, Real dt)
{
	


	//state_old_real.setDomainBndry(0, 0, 1, geom);
	//state_old_imag.setDomainBndry(0, 0, 1, geom);

    //FillDomainBoundary(state_old_real, geom, bc);
    //FillDomainBoundary(state_old_imag, geom, bc);

	getFunctional().evaluate(
		state_new_real,state_new_imag,
		state_old_real,state_old_imag,
		time);

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
	Real time, Real dt)
{
	

	evaluate(
		tmp_real,tmp_imag,
		state_old_real,state_old_imag,
		time);

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

	evaluate(
		tmp2_real,tmp2_imag,
		tmp_real,tmp_imag,
		time + dt*0.5
		);
	
	tmp2_real.mult(-dt);
	amrex::MultiFab::Saxpy(state_new_real,1./3.,tmp2_real,0,0,nComponents,ghosts);
	tmp2_imag.mult(-dt);
	amrex::MultiFab::Saxpy(state_new_imag,1./3.,tmp2_imag,0,0,nComponents,ghosts);

	tmp2_real.mult(0.5);
	tmp2_real.plus(state_old_real,0,nComponents,ghosts[0]);
	tmp2_imag.mult(0.5);
	tmp2_imag.plus(state_old_imag,0,nComponents,ghosts[0]);
	evaluate(
		tmp_real,tmp_imag,
		tmp2_real,tmp2_imag,
		time + dt*0.5);
	
	tmp_real.mult(-dt);
	tmp_imag.mult(-dt);
	amrex::MultiFab::Saxpy(state_new_real,1./3.,tmp_real,0,0,nComponents,ghosts);
	amrex::MultiFab::Saxpy(state_new_imag,1./3.,tmp_imag,0,0,nComponents,ghosts);

	tmp_real.plus(state_old_real,0,nComponents,ghosts[0]);
	tmp_imag.plus(state_old_imag,0,nComponents,ghosts[0]);
	evaluate(
		tmp2_real,tmp2_imag,
		tmp_real,tmp_imag,
		time + dt
		);
	tmp2_real.mult(-dt);
	tmp2_imag.mult(-dt);

	amrex::MultiFab::Saxpy(state_new_real,1./6.,tmp2_real,0,0,nComponents,ghosts);
	amrex::MultiFab::Saxpy(state_new_imag,1./6.,tmp2_imag,0,0,nComponents,ghosts);
}


void RK4Stepper::evaluate( 
	MultiFab & state_new_real, MultiFab & state_new_imag,
	MultiFab & state_old_real,  MultiFab & state_old_imag,
	Real time )

	
{

	auto & geom = getFunctional().getGeometry();

	state_old_real.setDomainBndry(0, 0, 1, geom);
	state_old_imag.setDomainBndry(0, 0, 1, geom);

    //FillDomainBoundary(state_old_real, geom, bc);
    //FillDomainBoundary(state_old_imag, geom, bc);

	
	getFunctional().evaluate(
			state_new_real,state_new_imag,
			state_old_real,state_old_imag,
			time
			);

	if (!isImaginaryTime)
	{
		state_new_imag.mult(-1);
		std::swap(state_new_real,state_new_imag);
	}

}

RK4Stepper::RK4Stepper(functional * func_  ,bool isImaginaryTime_,int nComponents_,int ghosts_) :
		ghosts(AMREX_D_DECL(ghosts_,ghosts_,ghosts_) ),nComponents(nComponents_),isImaginaryTime(isImaginaryTime_), stepper::stepper(func_) {

		auto & box = getFunctional().getBoxArray();
		auto & dm = getFunctional().getDistributionMapping();
		auto & geom = getFunctional().getGeometry();

		tmp_real.define(box, dm, nComponents, ghosts[0]);
		tmp2_real.define(box, dm , nComponents, ghosts[0]);
		tmp_imag.define( box, dm , nComponents, ghosts[0]);
		tmp2_imag.define( box, dm, nComponents, ghosts[0]);

		}

