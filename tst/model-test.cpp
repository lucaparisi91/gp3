#include "gtest/gtest.h"
#include "model.h"
#include <vector>
#include "evaluate.h"
#include "tools.h"
#include <cmath>
#include "operators.h"
#include "functional.h"
#include "run.h"
#include "functionalFactory.h"
#include "initializer.h"

TEST(initialization, createGeometry)
{
    
    auto settings = R"( 
        { 
        "geometry" : {"shape" : [128,128, 128] , "domain" : [ [-5,5] , [-5,5] , [-5,5] ] } } )"_json;
    
    auto [box , geom , dm] = createGeometry(settings["geometry"]);

    int Ncomp = 1;
    int Nghost = 2;

    MultiFab phi_real(box, dm, Ncomp, Nghost);

}


TEST(initialization, harmonic)
{
    auto settings = R"( 
        { 
        "geometry" : {"shape" : [128,128, 128] , "domain" : [ [-5,5] , [-5,5] , [-5,5] ] } } )"_json;
    
    auto [box , geom , dm] = createGeometry(settings["geometry"]);

    int Ncomp = 1;
    int Nghost = 2;


    MultiFab phi_real_old(box, dm, Ncomp, Nghost);
    MultiFab phi_imag_old(box, dm, Ncomp, Nghost);

    MultiFab phi_real_new(box, dm, Ncomp, Nghost);
    MultiFab phi_imag_new(box, dm, Ncomp, Nghost);

    phi_imag_old=0;

     auto functionalSettings = R"( 
        {
            "name" : "harmonic",
            "omega" : 1.0 
         } )"_json;
   
    auto func = initializer::instance().getFunctionalFactory().create(functionalSettings);

    func->define(geom,box,dm);

    Real alpha=1;
#if AMREX_SPACEDIM == 3
    fill(phi_real_old,geom, [alpha](Real x, Real y, Real z) {return exp(- alpha *(x*x + y*y + z*z));}  ) ;
#endif

#if AMREX_SPACEDIM == 1
    fill(phi_real_old,geom, [alpha](Real x) {return exp(- alpha *(x*x ));}  ) ;
#endif



    LOOP(phi_real_old,geom)

#if AMREX_SPACEDIM==3
     auto x = prob_lo[0] +  (i + 0.5) * dx[0] ;
	 auto y= prob_lo[1] + (j + 0.5) * dx[1];
	 auto z= prob_lo[2] + (k + 0.5) * dx[2];		

    auto r2 = x*x + y*y + z*z;

    ASSERT_NEAR(data(i,j,k), exp(-alpha*r2),1e-5 ) ;

#endif

#if AMREX_SPACEDIM==1
     auto r = prob_lo[0] +  (i + 0.5) * dx[0] ;
    ASSERT_NEAR(data(i,j,k), exp(-alpha*r*r),1e-5 ) ;
#endif

    ENDLOOP

    func->evaluate(phi_real_new,phi_imag_new,phi_real_old,phi_imag_old,0);

    LOOP(phi_real_new,geom)

#if AMREX_SPACEDIM == 3

     auto x = prob_lo[0] +  (i + 0.5) * dx[0] ;
	 auto y= prob_lo[1] + (j + 0.5) * dx[1];
	 auto z= prob_lo[2] + (k + 0.5) * dx[2];		

    auto r2 = x*x + y*y + z*z;

     ASSERT_NEAR( data(i,j,k,0) , exp(- alpha*r2 ) * ( 3*alpha  + r2 * (-2*alpha*alpha + 0.5 ) )   , 1e-2);
#endif

#if AMREX_SPACEDIM == 1

     auto x = prob_lo[0] +  (i + 0.5) * dx[0] ;
	
    auto r2 = x*x;

     ASSERT_NEAR( data(i,j,k,0) , exp(- alpha*r2 ) * ( alpha  + r2 * (-2*alpha*alpha + 0.5 ) )   , 1e-2);
#endif


    ENDLOOP


    delete func;

}