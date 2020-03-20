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
	virtual void evolve( 
	MultiFab & state_new_real, MultiFab & state_new_imag,
	MultiFab & state_old_real,  MultiFab & state_old_imag,
	MLPoisson & laplacianOperatorReal, MLPoisson & laplacianOperatorImag ,
	 Geometry & geom ,    const Vector<BCRec> & bc,
	 Real time, Real dt);

 };

 #endif