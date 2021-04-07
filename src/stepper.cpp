#include "stepper.h"
#include "timers.h"
namespace gp{

void eulerStepper::evolve( 
	wavefunction & waveNew, const wavefunction & waveOld,
	Real time, Real dt)
{
	//state_old_real.setDomainBndry(0, 0, 1, geom);
	//state_old_imag.setDomainBndry(0, 0, 1, geom);

    //FillDomainBoundary(state_old_real, geom, bc);
    //FillDomainBoundary(state_old_imag, geom, bc);


	auto & phi = waveNew.getPhi();
	const auto & phiOld = waveOld.getPhi();
	START_TIMER("evaluate");
	getFunctional().evaluate(
		waveNew, waveOld,
		time);
	STOP_TIMER("evaluate");


	if (!isImaginaryTime)
	{
		phi.mult(-dt);
		for(int c=0;c<waveNew.nComp()/2;c+=2 )
		{
			amrex::MultiFab::Swap (phi, phi, c, c+1, 1, phi.nGrow());

		}
	}
	else
	{
		phi.mult(-dt);
	}

	phi.plus(phiOld,0,phi.nComp(),0);
	START_TIMER("fill");
	waveNew.fillBoundaries();
	STOP_TIMER("fill");


}


void RK4Stepper::evolve( 
	wavefunction & waveNew, const wavefunction & waveOld
,
	Real time, Real dt)
{

	auto & phiNew = waveNew.getPhi();
	auto & phiTmp = waveTmp.getPhi();
	auto & phi2Tmp = waveTmp2.getPhi();

	const auto & phiOld = waveOld.getPhi();



	evaluate(
		waveTmp,waveOld,
		time);

	phiTmp.mult(-dt);

	
	amrex::MultiFab::Copy(phiNew, phiOld, 0, 0, nComponents,ghosts);
	amrex::MultiFab::Saxpy(phiNew,1./6.,phiTmp,0,0,nComponents,ghosts);


	phiTmp.mult(0.5);
	phiTmp.plus(phiOld,0,nComponents,ghosts[0]);


	evaluate(
	waveTmp2, waveTmp,
		time + dt*0.5
		);
	
	phi2Tmp.mult(-dt);
	amrex::MultiFab::Saxpy(phiNew,1./3.,phi2Tmp,0,0,nComponents,ghosts);
	

	phi2Tmp.mult(0.5);
	phi2Tmp.plus(phiOld,0,nComponents,ghosts[0]);
	
	evaluate(
		waveTmp, waveTmp2,
		time + dt*0.5);
	
	phiTmp.mult(-dt);
	

	amrex::MultiFab::Saxpy(phiNew,1./3.,phiTmp,0,0,nComponents,ghosts);
	
	phiTmp.plus(phiOld,0,nComponents,ghosts[0]);
	
	evaluate(
	waveTmp2, waveTmp,
		time + dt
		);
	phi2Tmp.mult(-dt);

	amrex::MultiFab::Saxpy(phiNew,1./6.,phi2Tmp,0,0,nComponents,ghosts);

}


void RK4Stepper::evaluate( 
	wavefunction & waveNew, const wavefunction & waveOld,
	Real time )	
{
	auto & phi = waveNew.getPhi();

	getFunctional().evaluate(
			waveNew, waveOld,
			time
			);

	if (!isImaginaryTime)
	{
		phi.mult(-1);

		for(int c=0;c<waveNew.nComp()/2;c+=2 )
		{
			amrex::MultiFab::Swap (phi, phi, c, c+1, 1, phi.nGrow());

		}
		
	}

	
	waveNew.fillBoundaries();

}

RK4Stepper::RK4Stepper(functional * func_  , bool isImaginaryTime_, wavefunction tmp_ , wavefunction tmp2_ ) :
stepper::stepper(func_),
waveTmp(tmp_),
waveTmp2(tmp2_),
isImaginaryTime(isImaginaryTime_)
  {

		}

}
