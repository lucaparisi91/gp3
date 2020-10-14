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
#include <AMReX_Orientation.H>

/*

TEST(twoLevel, laplacian)
{
    int nComp=1;
    int nGhost = 1;
     auto settingsCoarse = R"( 
        { 
        "geometry" : {"shape" : [128,128, 128] , "domain" : [ [-5,5] , [-5,5] , [-5,5] ] , "coordinates" : "cartesian"} } )"_json;

    auto settingsFine =  R"( 
        { 
        "geometry" : {"shape" : [256,256, 256] , "domain" : [ [-5,5] , [-5,5] , [-5,5] ] ,  "coordinates" : "cartesian" } } )"_json;


    auto [boxCoarse , geom , dmCoarse, low_bc, high_bc] = createGeometry(settingsCoarse["geometry"]);
    auto [boxFine , geomFine , dmFine, low_bcFine, high_bcFine] = createGeometry(settingsCoarse["geometry"]);

    Box innerDomain( {30,30,30},{60,60,60});
    BoxArray baInner;

    baInner.define(innerDomain);

    DistributionMapping dmInner(baInner);

    MultiFab coarseMultiFab( boxCoarse,dmCoarse,nComp,nGhost);
    MultiFab coarseTmp( boxCoarse,dmCoarse,nComp,nGhost);

    
    MultiFab innerMultiFab(baInner,dmInner,nComp,nGhost);

    innerMultiFab=0.;
    coarseTmp=0.;


    Real alpha = 0.01;

    fill( coarseMultiFab,geom, [alpha](Real x,Real y , Real z) {return exp(- alpha *(x*x + y*y + z*z ));}  ) ;

    // creates and fill a coarse boundary register
    auto baInnerCoarsened = baInner.coarsen(2);
    amrex::DistributionMapping dmInnerCoarsened(baInner);

    BndryRegister coarseBoundaryRegister(baInnerCoarsened,dmInnerCoarsened,0,1,0,1);
    coarseBoundaryRegister.copyFrom(coarseMultiFab, 1, 0, 0, 1 ,geom.periodicity());

    Orientation lowXOr( 0, amrex::Orientation::Side::low);

    coarseBoundaryRegister[lowXOr].copyTo(coarseTmp,1,0,0,1, geom.periodicity());

     LOOP(coarseTmp,geom)
    
    if ( data(i,j,k) != 0)
    {
        std::cout << data(i,j,k) << std::endl;
    }

    
     ENDLOOP






    // creates a fine boundry register
    
    auto bc =toMultiFabBC( low_bc,high_bc);
    amrex::InterpBndryData fineBoundary(baInner,dmInner,1,geomFine );
    fineBoundary.setBndryValues (coarseBoundaryRegister, 0,  innerMultiFab , 0, 0, 1, {2,2,2}, bc);
    

    fineBoundary[lowXOr].copyTo(innerMultiFab,1,0,0,1, geomFine.periodicity());

     LOOP( innerMultiFab,geom)
    
    if ( data(i,j,k) != 0)
    {
        std::cout << data(i,j,k) << std::endl;
    }
    
     ENDLOOP








    
    




}

*/