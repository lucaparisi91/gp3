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
#include <AMReX_Interpolater.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_Amr.H>
#include <AMReX_AmrLevel.H>
#include "AmrCoreDiff.h"



TEST(twoLevel, laplacian)
{
    Vector<int> is_periodic(AMREX_SPACEDIM,1);
    
    auto settingsCoarse = R"( 
        { 
        "geometry" : {"shape" : [20,20, 20] , "domain" : [ [-5,5] , [-5,5] , [-5,5] ] , "coordinates" : "cartesian"} } )"_json;

    auto settingsFine = R"( 
        { 
        "geometry" : {"shape" : [40,40, 40] , "domain" : [ [-5,5] , [-5,5] , [-5,5] ] , "coordinates" : "cartesian"} } )"_json;

    auto [baCoarse , geomCoarse , dmCoarse, low_bcCoarse, high_bcCoarse] = createGeometry(settingsCoarse["geometry"]);
    auto [baFine , geomFine , dmFine, low_bcFine, high_bcFine] = createGeometry(settingsFine["geometry"]);

    IntVect dom_lo(AMREX_D_DECL(       19 ,  19      ,  19  ));
    IntVect dom_hi(AMREX_D_DECL(       28,   28 ,  28));

    Box innerBox(dom_lo,dom_hi) ;
    BoxArray innerBoxArray;
    DistributionMapping dmInner(innerBoxArray);

    innerBoxArray.define(innerBox);
    int order = 2;
    int refinementRatio=2;
    int nComp=1;
    Vector<int> periodicity(AMREX_SPACEDIM,1);


    // fill the coarse multifab width a gaussian
    Real alpha=1.;
    MultiFab mfFine(baFine,dmFine,nComp,order-1);
    MultiFab mfCoarse(baCoarse,dmCoarse,nComp,order-1);


     for ( MFIter mfi( mfCoarse); mfi.isValid(); ++mfi )
    {
         const Real* dx = geomCoarse.CellSize();
        const Real* prob_lo = geomCoarse.ProbLo();

        const Box& bx = mfi.validbox();
        const int* lo = bx.loVect();
        const int *hi= bx.hiVect();
        Array4<  Real> const & dataArray = mfCoarse[mfi].array();
        
        for (int k=lo[2];k<=hi[2];k++)
            for (int j=lo[1];j<=hi[1];j++) 
                for (int i=lo[0];i<=hi[0];i++)
                    {
                        Real x =  prob_lo[0] + (i+0.5)*dx[0];
                        Real y =  prob_lo[1] + (i+0.5)*dx[1];
                        Real z =  prob_lo[2] + (i+0.5)*dx[2];

                        Real r2=x*x + y*y + z*z;
                        dataArray(i,j,k,0)=exp(-alpha*r2);
                    }
    }

    mfFine=0;

    auto bcRec = toMultiFabBC(low_bcCoarse,high_bcCoarse);

    amrex::CellQuadratic interp;

    //amrex::FillPatchTwoLevels(mfCoarse,{1,1,1},0.0,{&mfCoarse},{0.},{&mfFine},{0.0},0,0,1,geomCoarse,geomFine,BC::PERIODIC,0,BC::PERIODIC,0,{2,2,2},&interp,{bcRec},0);

    IntVect ghosts{1,1,1};
    
}


TEST(run, amrLevel)
{

    std::array<Real,AMREX_SPACEDIM> lower_edges={-3,-3,-3};
    std::array<Real,AMREX_SPACEDIM> higher_edges={3,3,3};


    RealBox real_box({AMREX_D_DECL( lower_edges[0],lower_edges[1],lower_edges[2]) },
                         {AMREX_D_DECL( higher_edges[0], higher_edges[1], higher_edges[2] )});
    
    Vector<int> nCells {32,32,32};

    int max_level = 1;

    AmrCoreDiff amr_core_diff(& real_box,max_level,nCells,0);

    amr_core_diff.InitData();

    amr_core_diff.ComputeDt();

    amr_core_diff.WritePlotFile();
    
    amr_core_diff.timeStepWithSubcycling(0, 0, 1);


    //amr_core_diff.AdvancePhiAtLevel(0,0,1e-3,0,1);
    //amr_core_diff.AdvancePhiAtLevel(1,0,1e-3,0,1);

  
    amr_core_diff.WritePlotFile();

                        
}