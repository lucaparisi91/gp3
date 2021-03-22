
#include "gpDropletFunctional3D.h"

void gpDroplet3D::evaluate(
	MultiFab & state_new_real, MultiFab & state_new_imag,
	MultiFab & state_old_real,  MultiFab & state_old_imag , Real time)
    {
        auto & lap = getLaplacian();

        lap.apply(state_new_real,state_old_real);
        lap.apply(state_new_imag,state_old_imag);

        auto & geom = getGeometry();

       EVALUATION_LOOP(state_new_real,state_new_imag,state_old_real,state_old_imag, getGeometry())


       phi_new_real(   i,j,k,0) *= -0.5 ;
	   phi_new_imag(  i,j,k,0) *= -0.5 ;
       
        
       Real psi2 = phi_old_real(i,j,k,0) * phi_old_real(i,j,k,0) + phi_old_imag(i,j,k,0) * phi_old_imag(i,j,k,0)  ;
       Real psi3 = std::pow(psi2,3/2.);

       Real intTerm = -3 * psi2 + 5/2. * psi3;

	   phi_new_real(  i,j,k,0)+= intTerm *phi_old_real(i,j,k,0);
	   phi_new_imag( i,j,k,0)+= intTerm*phi_old_imag( i,j,k,0);

       END_EVALUATION_LOOP
       
    } 