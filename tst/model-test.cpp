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

    functionalFactory facFunc;
    
    facFunc.registerFunctional<harmonicFunctional>();

     auto functionalSettings = R"( 
        {
            "omega" : 1.0 
         } )"_json;
   
    auto func = facFunc.create("harmonic",functionalSettings);

    func->define(geom,box,dm);

    Real alpha=0.87;

    fill(phi_real_old,geom, [alpha](Real x, Real y, Real z) {return exp(- alpha *(x*x + y*y + z*z));}  ) ;

    LOOP3D(phi_real_old,geom)

     auto x = prob_lo[0] +  (i + 0.5) * dx[0] ;
	 auto y= prob_lo[1] + (j + 0.5) * dx[1];
	 auto z= prob_lo[2] + (k + 0.5) * dx[2];		

    auto r2 = x*x + y*y + z*z;

    ASSERT_NEAR(data(i,j,k), exp(-alpha*r2),1e-5 ) ;


    ENDLOOP3D

    func->evaluate(phi_real_new,phi_imag_new,phi_real_old,phi_imag_old,0);

    LOOP3D(phi_real_new,geom)

     auto x = prob_lo[0] +  (i + 0.5) * dx[0] ;
	 auto y= prob_lo[1] + (j + 0.5) * dx[1];
	 auto z= prob_lo[2] + (k + 0.5) * dx[2];		

    auto r2 = x*x + y*y + z*z;

     ASSERT_NEAR( data(i,j,k,0) , exp(- alpha*r2 ) * ( 3*alpha  + r2 * (-2*alpha*alpha + 0.5 ) )   , 1e-2);

    ENDLOOP3D



    delete func;

}