#ifndef STEPPER_H
#define STEPPER_H

#include <AMReX_MultiFabUtil.H>
#include <AMReX_MLMG.H>
#include "evaluate.h"
#include "functional.h"

using namespace amrex;

class stepper
 {

	 
public:

	stepper( functional * func_) : _func(func_) {}

	virtual void evolve( 
	MultiFab & state_new_real, MultiFab & state_new_imag,
	MultiFab & state_old_real,  MultiFab & state_old_imag,
	 Real time, Real dt)=0;

	auto & getFunctional() {return *_func;} 

private:

	 functional * _func;
 };


 class eulerStepper : public stepper
 {
 public:
 	eulerStepper(functional * func_,bool isImaginaryTime_) : isImaginaryTime(isImaginaryTime_) , stepper::stepper(func_){}


	virtual void evolve( 
	MultiFab & state_new_real, MultiFab & state_new_imag,
	MultiFab & state_old_real,  MultiFab & state_old_imag,
	 Real time, Real dt);
 private:
 	bool isImaginaryTime;

 };


 class RK4Stepper : public stepper
 {
 	/* Implements the 4th order Runje-Kutta stepper*/
 public:
 	RK4Stepper(functional * func_  ,bool isImaginaryTime_,int nComponents_,int ghosts_);
	 
 	virtual void evolve( 
 	MultiFab & state_new_real, MultiFab & state_new_imag,
 	MultiFab & state_old_real,  MultiFab & state_old_imag,
 	 Real time, Real dt);
private:
	void evaluate( 
			MultiFab & state_new_real, MultiFab & state_new_imag,
			MultiFab & state_old_real,  MultiFab & state_old_imag,
			Real time );

	MultiFab tmp_real;
	MultiFab tmp_imag;
	MultiFab tmp2_real;
	MultiFab tmp2_imag;
	int nComponents;
	amrex::IntVect ghosts;
	
	bool isImaginaryTime;
 };

 #endif