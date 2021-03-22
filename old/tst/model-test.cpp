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
#include "plotFile.h"


TEST(initialization, createGeometry)
{
    
    auto settings = R"( 
        { 
        "geometry" : {"shape" : [128,128, 128] , "domain" : [ [-5,5] , [-5,5] , [-5,5] ] , "coordinates" : "cartesian"} } )"_json;
    
    auto [box , geom , dm, low_bc, high_bc] = createGeometry(settings["geometry"]);

    int Ncomp = 1;
    int Nghost = 2;

    MultiFab phi_real(box, dm, Ncomp, Nghost);
}

TEST(initialization, harmonic)
{
    auto settings = R"( 
        { 
        "geometry" : {"shape" : [128,128, 128] , "domain" : [ [-5,5] , [-5,5] , [-5,5] ],"coordinates" : "cartesian" } } )"_json;
    
    auto [box , geom , dm, low_bc, high_bc] = createGeometry(settings["geometry"]);

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
            "omega" : 1.0 ,
            "laplacian" : {
                "name" : "amrexLaplacian",
                "order" : 2
            }
         } )"_json;
    
    auto func = initializer::instance().getFunctionalFactory().create(functionalSettings);

    func->define(geom,box,dm,low_bc,high_bc);

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

     //ASSERT_NEAR( data(i,j,k,0) , exp(- alpha*r2 ) * ( alpha  + r2 * (-2*alpha*alpha + 0.5 ) )   , 1e-2);
#endif


    ENDLOOP


    delete func;

}


#if AMREX_SPACEDIM == 1

TEST(initialization,harmonic_spherical)
{
    auto settings = R"( 
        { 
        "geometry" : {"shape" : [1000] , "domain" : [ [0,10]  ],"coordinates" : "spherical" } } )"_json;
    
    auto [box , geom , dm, low_bc, high_bc] = createGeometry(settings["geometry"]);

    int Ncomp = 1;
    int Nghost = 1;


    MultiFab phi_real_old(box, dm, Ncomp, Nghost);
    MultiFab phi_imag_old(box, dm, Ncomp, Nghost);
    
    MultiFab phi_real_new(box, dm, Ncomp, Nghost);
    MultiFab phi_imag_new(box, dm, Ncomp, Nghost);


    Real alpha = 1;
    phi_real_old=0.;
    phi_real_new=0.;

    fill(phi_real_old,geom, [alpha](Real x) {return exp(- alpha *(x*x ));}  ) ;


      auto functionalSettings = R"( 
        {
            "name" : "harmonic",
            "omega" : 1.0 ,
            "laplacian" :
            {
                "name" : "stencilLaplacian",
                "order" : 2 
            }
         } )"_json;
    
    auto func = initializer::instance().getFunctionalFactory().create(functionalSettings);
    func->define(geom,box,dm,low_bc,high_bc);

    func->evaluate(phi_real_new,phi_imag_new,phi_real_old,phi_imag_old,0);

    LOOP(phi_real_new,geom)

     auto r = prob_lo[0] +  (i + 0.5) * dx[0] ;
	 
     auto r2 = r*r;


     ASSERT_NEAR( data(i,j,k,0) , exp(- alpha*r2 ) * ( 3*alpha  + r2 * (-2*alpha*alpha + 0.5 ) )   , 1e-2);

    ENDLOOP


    auto norm2 = norm(phi_real_old,phi_imag_old,geom);

    ASSERT_NEAR(norm2 , std::pow(M_PI/(2*alpha), 3/2.) ,1e-2);

    normalize(phi_real_old,phi_imag_old , geom);

    norm2 = norm(phi_real_old,phi_imag_old,geom);

    ASSERT_NEAR(norm2,1,1e-2);


}


#endif


TEST(initialization, harmonic_stencil2)
{
    
    auto settings = R"( 
        { 
        "geometry" : {"shape" : [200 , 200 , 200] , "domain" : [ [-5,5] , [-5,5] , [-5,5] ],"coordinates" : "cartesian" } } )"_json;
    
    auto [box , geom , dm, low_bc, high_bc] = createGeometry(settings["geometry"]);

    int Ncomp = 1;
    int Nghost = 1;

    MultiFab phi_real_old(box, dm, Ncomp, Nghost);
    MultiFab phi_imag_old(box, dm, Ncomp, Nghost);

    MultiFab phi_real_new(box, dm, Ncomp, Nghost);
    MultiFab phi_imag_new(box, dm, Ncomp, Nghost);

    phi_imag_old=0;
    phi_real_old=0;

     auto functionalSettings = R"( 
        {
            "name" : "harmonic",
            "omega" : 1.0 ,
            "laplacian" : {
                "name" : "stencilLaplacian2",
                "order" : 2
            }
         } )"_json;
    
    auto func = initializer::instance().getFunctionalFactory().create(functionalSettings);

    func->define(geom,box,dm,low_bc,high_bc);

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



#if AMREX_SPACE_DIM == 1

TEST(initialization, harmonic_stencil2_spherical)
{
    auto settings = R"( 
        { 
        "geometry" : {"shape" : [1000] , 
        "domain" : [ [0,5] ],
        "coordinates" : "spherical" ,
        "bc" : ["neumann"]
        } } )"_json;
        
    auto [box , geom , dm, low_bc, high_bc] = createGeometry(settings["geometry"]);

    int Ncomp = 1;
    int Nghost = 1;

    MultiFab phi_real_old(box, dm, Ncomp, Nghost);
    MultiFab phi_imag_old(box, dm, Ncomp, Nghost);

    MultiFab phi_real_new(box, dm, Ncomp, Nghost);
    MultiFab phi_imag_new(box, dm, Ncomp, Nghost);

    phi_imag_old=0;
    phi_real_old=0;

     auto functionalSettings = R"( 
        {
            "name" : "harmonic",
            "omega" : 1.0 ,
            "laplacian" : {
                "name" : "stencilLaplacian2",
                "order" : 2
            }
         } )"_json;
    
    auto func = initializer::instance().getFunctionalFactory().create(functionalSettings);

    func->define(geom,box,dm,low_bc,high_bc);
    Real alpha = 1.;
    fill(phi_real_old,geom, [alpha](Real x) {return exp(- alpha *(x*x ));}  ) ;

    phi_real_old.setDomainBndry(0, 0, 1, geom)   ;
    phi_imag_old.setDomainBndry(0, 0, 1, geom)   ;


    LOOP(phi_real_old,geom)

     auto r = prob_lo[0] +  (i + 0.5) * dx[0] ;
    ASSERT_NEAR(data(i,j,k), exp(-alpha*r*r),1e-5 ) ;


    ENDLOOP

    func->evaluate(phi_real_new,phi_imag_new,phi_real_old,phi_imag_old,0);

    LOOP(phi_real_new,geom)


     auto x = prob_lo[0] +  (i + 0.5) * dx[0] ;

    auto r2 = x*x;

    {
    ASSERT_NEAR( data(i,j,k,0) , exp(- alpha*r2 ) * ( 3* alpha  + r2 * (-2*alpha*alpha + 0.5 ) )   , 1e-2);
    }

    ENDLOOP


    delete func;

}

#endif



TEST(save, parquetTest)
{
    #if AMREX_SPACEDIM == 1

    auto settings = R"( 
        { 
        "geometry" : {"shape" : [1000] , 
        "domain" : [ [0,5] ],
        "coordinates" : "spherical" ,
        "bc" : ["neumann"]
        } } )"_json;
    #endif


    #if AMREX_SPACEDIM == 3

    auto settings = R"( 
        { 
        "geometry" : {"shape" : [100,100,100] , 
        "domain" : [ [-5,5],[-5,5],[-5,5] ],
        "coordinates" : "cartesian" ,
        "bc" : ["periodic","periodic","periodic"]
        } } )"_json;

    #endif


        
    auto [box , geom , dm, low_bc, high_bc] = createGeometry(settings["geometry"]);

    int Ncomp = 1;
    int Nghost = 1;

    MultiFab phi_real_old(box, dm, Ncomp, Nghost);
    MultiFab phi_imag_old(box, dm, Ncomp, Nghost);

    MultiFab phi_real_new(box, dm, Ncomp, Nghost);
    MultiFab phi_imag_new(box, dm, Ncomp, Nghost);

    phi_imag_old=0;
    phi_real_old=0;

    Real alpha = 1;

    #if AMREX_SPACEDIM == 1
    fill(phi_real_old,geom, [alpha](Real x) {return exp(- alpha *(x*x ));}  ) ;
    #endif

    #if AMREX_SPACEDIM == 3
    fill(phi_real_old,geom, [alpha](Real x,Real y , Real z) {return exp(- alpha *(x*x + y*y + z*z ));}  ) ;
    #endif




    phi_real_old.setDomainBndry(0, 0, 1, geom)   ;
    phi_imag_old.setDomainBndry(0, 0, 1, geom)   ;

    auto bc = toMultiFabBC(low_bc,high_bc) ;

    FillDomainBoundary(phi_real_old, geom, {bc} );



    LOOP(phi_real_old,geom)
    #if AMREX_SPACE_DIM == 1
     auto r = prob_lo[0] +  (i + 0.5) * dx[0] ;
     #endif

    #if AMREX_SPACEDIM == 3
     auto x = prob_lo[0] +  (i + 0.5) * dx[0] ;
     auto y = prob_lo[0] +  (j + 0.5) * dx[1] ;
     auto z = prob_lo[0] +  (k + 0.5) * dx[2] ;

     auto r = std::sqrt( x*x + y*y + z*z);
     #endif
     
    ASSERT_NEAR(data(i,j,k), exp(-alpha*r*r),1e-5 ) ;
    ENDLOOP

    writeSingleLevel(phi_real_old,phi_imag_old,geom, "out-test");

}