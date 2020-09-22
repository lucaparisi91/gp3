#include "gtest/gtest.h"
#include "model.h"
#include <vector>
#include "evaluate.h"
#include "tools.h"
#include <cmath>

TEST(initModel,evaluate)
{

    std::array<size_t,AMREX_SPACEDIM> shape { 300 , 300 , 300};
    std::array<Real,AMREX_SPACEDIM> lower_edges {-5 , -5 ,-5 };
    std::array<Real,AMREX_SPACEDIM> higher_edges {5 , 5 ,5 };

    geometry geom( shape );

    geom.lower_edges = lower_edges;
    geom.higher_edges = higher_edges;

    int order = 4;



    model m(geom,order, 1);
    model m2(geom,order, 1);

    
    Real alpha =0.5;
    

    m.fill( [alpha](Real x, Real y, Real z) {return exp(- alpha *(x*x + y*y + z*z));}  );

    auto initNorm2 = norm( m.real() , m.imag(),  m.getGeometry(), 0) ;

    Real initNorm2Expected = std::pow( (M_PI/(2*alpha)) , 3/2. );

    ASSERT_NEAR(initNorm2 , initNorm2Expected, 1e-2 );

    auto & state = m.real();

    LOOP( state , m.getGeometry() )

    auto x = prob_lo[0] +  (i + 0.5) * dx[0] ;
	auto y= prob_lo[1] + (j + 0.5) * dx[1];
	auto z= prob_lo[2] + (k + 0.5) * dx[2];		

    ASSERT_NEAR( data(i,j,k,0) , exp(-alpha* ( x*x +  y*y + z*z) )  , 1e-2);

    ENDLOOP


    evaluate(m2.real(), m2.imag() ,
        m.real() , m.imag() , 
        0 , m.getGeometry() , m.realLaplacian() , m.imagLaplacian() 
     );




    auto & state2 = m2.real();
    
    LOOP( state2 , m2.getGeometry() )

        auto x = prob_lo[0] +  (i + 0.5) * dx[0] ;
	    auto y= prob_lo[1] + (j + 0.5) * dx[1];
	    auto z= prob_lo[2] + (k + 0.5) * dx[2];		

        auto r2 = x*x + y*y + z*z;
         ASSERT_NEAR( data(i,j,k,0) , exp(- alpha*r2 ) * ( 3*alpha  + r2 * (-2*alpha*alpha + 0.5) )   , 1e-2);
         //ASSERT_NEAR( data(i,j,k,0) , exp(-alpha* ( x*x +  y*y + z*z) )  , 1e-2);

         //ASSERT_NEAR( data(i,j,k,0) , 0  , 1e-2);


    ENDLOOP







}