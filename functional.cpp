#include "functional.h"

void functional::define( Geometry & geom_ , BoxArray & box_, DistributionMapping & dm_ )
{
    _geom = geom_;
    _box = box_;
    _dm=dm_;

    // creates a new laplacian if one was not previously set
    if (lap == NULL)
    {
        lap=new amrexLaplacianOperator();
        laplacianOwned=true;
    }

    lap->define(geom_,box_,dm_);


};

functional::~functional()
{
    if (laplacianOwned)
    {delete lap;}
}

void functional::setLaplacianOperator(laplacianOperator * lap_)
{
    if (laplacianOwned) 
    {
        delete lap;
    };

    lap=lap_;
    laplacianOwned=false;
}

 void harmonicFunctional::evaluate(
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


	    auto [ AMREX_D_DECL( x , y , z) ] = std::tuple(AMREX_D_DECL( prob_lo[0] +  (i + 0.5) * dx[0] ,
	         prob_lo[1] + (j + 0.5) * dx[1] , 
	   	    prob_lo[2] + (k + 0.5) * dx[2] ) );
        
        #if AMREX_SPACEDIM == 3
        Real tmp =   0.5 *( x*x + y*y + z*z);
        #endif

        #if AMREX_SPACEDIM == 1
        Real tmp =   0.5 *( x*x);
        #endif
    
		phi_new_real(  i,j,k,0)+= tmp*phi_old_real(i,j,k,0);
		phi_new_imag( i,j,k,0)+= tmp*phi_old_imag( i,j,k,0);


       END_EVALUATION_LOOP




    } 