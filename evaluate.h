#include <AMReX_BCUtil.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_MLPoisson.H>

using namespace amrex;
/*
Evaluates the new state given the old state and the current time 
and the laplacian operator
*/
//void evaluate( MultiFab & data_new,  const MultiFab & data_old, Real time,  MLMG & laplacianOperator );

void evaluate( 
	MultiFab & state_new_real, MultiFab & state_new_imag,
	MultiFab & state_old_real,  MultiFab & state_old_imag,
	Real time, Geometry & geom ,  MLPoisson & laplacianOperatorReal, MLPoisson & laplacianOperatorImag );