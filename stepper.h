#ifndef STEPPER_H
#define STEPPER_H

#include <AMReX_MultiFabUtil.H>
#include <AMReX_MLMG.H>
#include "evaluate.h"
using namespace amrex;


 class stepper
 {
public:
	virtual void evolve( 
	MultiFab & state_new_real, MultiFab & state_new_imag,
	MultiFab & state_old_real,  MultiFab & state_old_imag,
	MLPoisson & laplacianOperatorReal, MLPoisson & laplacianOperatorImag ,
	 Geometry & geom ,    const Vector<BCRec> & bc,
	 Real time, Real dt)=0;
 };

 class eulerStepper : public stepper
 {
 public:

 	eulerStepper(bool isImaginaryTime_) : isImaginaryTime(isImaginaryTime_) {}
	virtual void evolve( 
	MultiFab & state_new_real, MultiFab & state_new_imag,
	MultiFab & state_old_real,  MultiFab & state_old_imag,
	MLPoisson & laplacianOperatorReal, MLPoisson & laplacianOperatorImag ,
	 Geometry & geom ,    const Vector<BCRec> & bc,
	 Real time, Real dt);
 private:
 	bool isImaginaryTime;
 };

 class RK4Stepper : public stepper
 {
 	/* Implements the 4th order Runje-Kutta stepper*/
 public:
 	RK4Stepper(bool isImaginaryTime_,int nComponents_,int ghosts_) : 
 		ghosts({ghosts_,ghosts_,ghosts_}),nComponents(nComponents_),isImaginaryTime(isImaginaryTime_) {}
 	virtual void evolve( 
 	MultiFab & state_new_real, MultiFab & state_new_imag,
 	MultiFab & state_old_real,  MultiFab & state_old_imag,
 	MLPoisson & laplacianOperatorReal, MLPoisson & laplacianOperatorImag ,
 	 Geometry & geom ,    const Vector<BCRec> & bc,
 	 Real time, Real dt);
private:
	MultiFab tmp_real;
	MultiFab tmp_imag;
	MultiFab tmp2_real;
	MultiFab tmp2_imag;
	int nComponents;
	IntVect ghosts;
	void evaluate_complex_( 
	MultiFab & state_new_real, MultiFab & state_new_imag,
	MultiFab & state_old_real,  MultiFab & state_old_imag,
	Real time, Geometry & geom ,  MLPoisson & laplacianOperatorReal, MLPoisson & laplacianOperatorImag );
	bool isImaginaryTime;
 };

 #endif