#ifndef STEPPER_H
#define STEPPER_H

#include <AMReX_MultiFabUtil.H>
#include <AMReX_MLMG.H>
#include "functional.h"
#include "wavefunction.h"



using namespace amrex;

namespace gp
{

class stepper
 {

public:
	stepper( functional * func_) : _func(func_) {}

	virtual void evolve( 
		wavefunction & waveNew, const wavefunction & waveOld
	,
	 Real time, Real dt)=0;

	auto & getFunctional() {return *_func;} 

private:

	 functional * _func;
 };


 class eulerStepper : public stepper
 {
 public:
 	eulerStepper(functional * func_,bool isImaginaryTime_) : stepper::stepper(func_),isImaginaryTime(isImaginaryTime_) {}
	

	virtual void evolve( 
	wavefunction & waveNew, const wavefunction & waveOld,
	 Real time, Real dt);


 private:
 	bool isImaginaryTime;
 };

 class RK4Stepper : public stepper
 {
 	/* Implements the 4th order Runge-Kutta stepper*/
 public:
 	RK4Stepper(functional * func_  ,bool isImaginaryTime_,wavefunction  tmp, wavefunction tmp2);

	 
 	virtual void evolve( 
		 wavefunction & waveNew, const wavefunction & waveOld,
 	 Real time, Real dt);
private:

	void evaluate( 
			wavefunction & waveNew, const wavefunction & waveOld,
			Real time );


	wavefunction waveTmp;
	wavefunction waveTmp2;

	int nComponents;
	amrex::IntVect ghosts;

	bool isImaginaryTime;
 };

}
 #endif