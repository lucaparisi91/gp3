#include "functional.h"

void functional::define( Geometry & geom_ , BoxArray & box_, DistributionMapping & dm_ )
{
    _geom = geom_;
    _box = box_;
    _dm=dm_;

    lap=new amrexLaplacianOperator();
    lap->define(geom_,box_,dm_);



};

functional::~functional()
{
    if (lap != NULL)
    {delete lap;}
}

 void harmonicFunctional::evaluate(
	MultiFab & state_new_real, MultiFab & state_new_imag,
	MultiFab & state_old_real,  MultiFab & state_old_imag )
    {
        auto & lap = getLaplacian();

        lap.apply(state_new_real,state_old_real);
        lap.apply(state_new_imag,state_old_imag);

        auto & geom = getGeometry();

       EVALUATION_LOOP3D(state_new_real,state_new_imag,state_old_real,state_old_imag, getGeometry())

       phi_new_real(i,j,k,0) *= -0.5 ;
	   phi_new_imag(i,j,k,0) *= -0.5 ;

	    auto x = prob_lo[0] +  (i + 0.5) * dx[0] ;
	    auto y= prob_lo[1] + (j + 0.5) * dx[1];
	   	auto z= prob_lo[2] + (k + 0.5) * dx[2];
        Real tmp =   0.5 *( x*x + y*y + z*z);
		phi_new_real(i,j,k,0)+= tmp*phi_old_real(i,j,k,0);
		phi_new_imag(i,j,k,0)+= tmp*phi_old_imag(i,j,k,0);

       END_EVALUATION_LOOP3D




    } 