#include "evaluate.h"
#include "stencils.h"

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

	 	for (int j=lo[1];j<=hi[1];j++)
	     	for (int i=lo[0];i<=hi[0];i++)
	     	{
	     		phi_new_real(i,j,0) *= -0.5 ;
	    // 		 /*
	    // 	Non linear interactions
	    // 		*/	
	     		auto v = phi_old_real(i,j,0)*phi_old_real(i,j,0) +  
	     		         phi_old_imag(i,j,0)*phi_old_imag(i,j,0);
	     		phi_new_real(i,j,0) += 1. * v *phi_old_real(i,j,0) ;
	     		phi_new_real(i,j,0) += 1. * v * phi_old_imag(i,j,0) ;

	     		auto r = prob_lo[0] +  (i + 0.5) * dx[0] ;
	     		auto z= prob_lo[1] + (j + 0.5) * dx[1]; 
	    // 		/*
	    // 			Trapping potential
	    // 		*/	
	    	
	     		phi_new_real(i,j,0)+=0.5  *( 1. * r*r + 1. *z*z) * phi_old_real(i,j,0);
	     		phi_new_imag(i,j,0)+=0.5* (1. + r*r + 1. * z*z) * phi_old_imag(i,j,0);

	 	    }	

	 	

	 	
	}

	

}