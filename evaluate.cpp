#include "evaluate.h"
#include "stencils.h"
#include "tools.h"
#include <pybind11/pybind11.h>


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
	     	{ 

#define ENDLOOP \
			 }\
	} \
	}


void evaluate( 
	MultiFab & state_new_real, MultiFab & state_new_imag,
	MultiFab & state_old_real,  MultiFab & state_old_imag,
	Real time, Geometry & geom ,MLPoisson & laplacianOperatorReal ,MLPoisson & laplacianOperatorImag )
{
	// laplacian terms

    MLMG lapReal(laplacianOperatorReal);
    MLMG lapImag(laplacianOperatorImag);

	
	lapReal.apply({&state_new_real},{ &state_old_real});
	lapImag.apply({&state_new_imag},{ &state_old_imag});

	
	//state_new_real.mult(-0.5, 0, 0);
	//state_new_imag.mult(-0.5, 0, 0);

	const Real* dx = geom.CellSize();
	const Real* prob_lo = geom.ProbLo();

	for ( MFIter mfi(state_new_real); mfi.isValid(); ++mfi )
	{
	    const Box& bx = mfi.validbox();
	    const int* lo = bx.loVect();
	    const int *hi= bx.hiVect();
	    Array4< const Real> const & phi_old_real = state_old_real[mfi].const_array();
	    Array4< Real> const & phi_new_real = state_new_real[mfi].array();

	    Array4< const Real> const & phi_old_imag = state_old_imag[mfi].const_array();
	    Array4< Real> const & phi_new_imag = state_new_imag[mfi].array();

#if AMREX_SPACEDIM == 3
		for (int k=lo[2];k<=hi[2];k++)
	 		for (int j=lo[1];j<=hi[1];j++)
	     		for (int i=lo[0];i<=hi[0];i++)
	     	{
	     		phi_new_real(i,j,k,0) *= -0.5 ;
	     		phi_new_imag(i,j,k,0) *= -0.5 ;

	     		auto x = prob_lo[0] +  (i + 0.5) * dx[0] ;
	     		auto y= prob_lo[1] + (j + 0.5) * dx[1];
				auto z= prob_lo[2] + (k + 0.5) * dx[2];				
				// * (phi_old_real(i,j,k,0)*phi_old_real(i,j,k,0) + phi_old_imag(i,j,k,0)*phi_old_imag(i,j,k,0)) 
	     		Real tmp =   0.5 *( x*x + y*y + z*z);
				phi_new_real(i,j,k,0)+= tmp*phi_old_real(i,j,k,0);
				phi_new_imag(i,j,k,0)+= tmp*phi_old_imag(i,j,k,0);

	 	    }	
#endif
	 	

	 	
	}

	

}




void evaluate(model & modelNew,  model & modelOld, Real time)
{
	evaluate( 
	modelNew.real(), modelNew.imag(),
	modelOld.real(), modelOld.imag(),
	time, modelOld.getGeometry() ,modelOld.realLaplacian() ,
	modelOld.imagLaplacian()
	);
}