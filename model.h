#ifndef MODEL_H
#define MODEL_H

#include "geometry.h"
#include<array>
#include "traits.h"
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

using namespace amrex;

/* Fills the multifab with values from a 3D numpy array*/
void fillMultiFab( MultiFab & fab  , const py::array_t<double> & values );


/* Defines the model. 
Setub the box,distribution mapping 
Also builds the laplacian and gradient operators.
Defined only for a single level.
 */

class model
{
public:
    model(const geometry & geomInfo, int laplacianOrder=2, int NComp=1);
    
    MultiFab & real() {return phi_real;}
    MultiFab & imag() {return phi_imag;}

    const MultiFab & real() const {return phi_real;}

    auto & realLaplacian() {return linPoissonReal;}
    auto & imagLaplacian() {return linPoissonImag;}

    auto & getGeometry() {return geom;}


    void fill( const py::array_t<double> & real, const py::array_t<double> & imag  );

    template<class callable_t>
    void fill( const callable_t & evaluator ) 
    {

	const Real* dx = geom.CellSize();
	const Real* prob_lo = geom.ProbLo();

	for ( MFIter mfi(phi_real); mfi.isValid(); ++mfi )
	{
	    const Box& bx = mfi.validbox();
	    const int* lo = bx.loVect();
	    const int *hi= bx.hiVect();
	    Array4< Real> const & dataReal = phi_real[mfi].array();

#if AMREX_SPACEDIM == 3
		for (int k=lo[2];k<=hi[2];k++)
	 		for (int j=lo[1];j<=hi[1];j++)
	     		for (int i=lo[0];i<=hi[0];i++)
	     	{
	    
	     		auto x = prob_lo[0] +  (i + 0.5) * dx[0] ;
	     		auto y= prob_lo[1] + (j + 0.5) * dx[1];
				auto z= prob_lo[2] + (k + 0.5) * dx[2];

                dataReal( i, j , k , 0 )=  exp(-x*x - y*y - z*z);
                			
		
	 	    }	
#endif
	 	
}
}




private:
    BoxArray ba;
    Geometry geom;

    int laplacianOrder;
    int Nghost;
    int Ncomp;
    DistributionMapping dmInit;
    MultiFab phi_real;
    MultiFab phi_imag;

    MLPoisson linPoissonReal;
    MLPoisson linPoissonImag;

};




#endif