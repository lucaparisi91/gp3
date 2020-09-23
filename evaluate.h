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


#include "model.h"

#define LOOP( state , geom  ) \
	{ \
	const Real* dx = geom.CellSize(); \
	const Real* prob_lo = geom.ProbLo(); \
	for ( MFIter mfi(state); mfi.isValid(); ++mfi ) \
	{ \
	    const Box& bx = mfi.validbox(); \
	    const int* lo = bx.loVect(); \
	    const int *hi= bx.hiVect(); \
	    Array4< Real> const & data = state[mfi].array();\
		for (int k=lo[2];k<=hi[2];k++) \
	 		for (int j=lo[1];j<=hi[1];j++) \
	     		for (int i=lo[0];i<=hi[0];i++) \
	     	{ \

#define ENDLOOP \
			 }\
	} \
	} 

