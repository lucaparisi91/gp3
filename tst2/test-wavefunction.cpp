#include "gtest/gtest.h"
#include "geometry.h"
#include "wavefunction.h"
#include "kernelLauncher.h"
#include "functional.h"
#include "stepper.h"
#include "normalization.h"


struct check_t{
    gp::Real tol = 1e-5;
};


struct lapCheck_t {gp::Real tol = 1e-4; };

struct gpCheck_t {gp::Real tol = 5e-2; gp::Real g=0; gp::Real omega=1; };

struct gaussian : public gp::operators::operatorBase
{
    struct checkNormalization_t{
    gp::Real tol = 5e-2;
    };


    gaussian(const gp::Real alpha_) : alpha(alpha_)
    {
        
    };

    gp::Real operator()(int i, int j , int k, int c, gp::wavefunctionRegion & wave)
    {
       
        const auto & geom = wave.getGeometry();

        gp::Real x = geom.ProbLo()[0] + geom.CellSize()[0]*i  ;
        gp::Real y = geom.ProbLo()[1] + geom.CellSize()[1]*j  ;
        gp::Real z = geom.ProbLo()[2] + geom.CellSize()[2]*k  ;

        return exp(-alpha*(x*x + y*y + z*z)  );
           
    }


    void operator()(int i, int j , int k, int c, gp::wavefunctionRegion & wave,  const check_t & check)
    {
        auto & phi = wave.getPhi<gp::array_t>();

        const auto & geom = wave.getGeometry();

        gp::Real x = geom.ProbLo()[0] + geom.CellSize()[0]*i  ;
        gp::Real y = geom.ProbLo()[1] + geom.CellSize()[1]*j  ;
        gp::Real z = geom.ProbLo()[2] + geom.CellSize()[2]*k  ;

        ASSERT_NEAR( phi(i,j,k,c)  , exp(-alpha*(x*x + y*y + z*z))  , check.tol );

    }

    void operator()(int i, int j , int k, int c, gp::wavefunctionRegion & wave, const lapCheck_t & check)
    {
        auto & phi = wave.getPhi<gp::array_t>();

        const auto & geom = wave.getGeometry();

        gp::Real x = geom.ProbLo()[0] + geom.CellSize()[0]*i  ;
        gp::Real y = geom.ProbLo()[1] + geom.CellSize()[1]*j  ;
        gp::Real z = geom.ProbLo()[2] + geom.CellSize()[2]*k  ;
        auto r2 = x*x + y*y + z*z;


        ASSERT_NEAR( phi(i,j,k,c)  , exp(-alpha*r2) *(-2*alpha*(3-2*alpha*r2) )  , check.tol );

    }

    void operator()(int i, int j , int k, int c, gp::wavefunctionRegion & wave, const gpCheck_t & check)
    {
        auto & phi = wave.getPhi<gp::array_t>();

        const auto & geom = wave.getGeometry();

        gp::Real x = geom.ProbLo()[0] + geom.CellSize()[0]*i  ;
        gp::Real y = geom.ProbLo()[1] + geom.CellSize()[1]*j  ;
        gp::Real z = geom.ProbLo()[2] + geom.CellSize()[2]*k  ;
        auto r2 = x*x + y*y + z*z;

        ASSERT_NEAR( phi(i,j,k,c)  , -0.5*exp(-alpha*r2) *(-2*alpha*(3-2*alpha*r2) ) + 0.5*check.omega * r2 , check.tol );

    }

    void operator()(int i, int j , int k, int c , gp::wavefunctionRegion & wave, const checkNormalization_t & check)
    {
        auto & phi = wave.getPhi<gp::array_t>();
        const auto & geom = wave.getGeometry();


        gp::Real x = geom.ProbLo()[0] + geom.CellSize()[0]*i  ;
        gp::Real y = geom.ProbLo()[1] + geom.CellSize()[1]*j  ;
        gp::Real z = geom.ProbLo()[2] + geom.CellSize()[2]*k  ;

        gp::Real sigma2 = 1.0/(4*alpha);

        ASSERT_NEAR( phi(i,j,k,c)  , exp(-alpha*(x*x + y*y + z*z))/( std::pow(2*M_PI*sigma2,3/4.) * std::sqrt(2) )  , check.tol );

    }


    private:
    const gp::Real alpha;

};





TEST(wavefunction, evaluateFunctional)
{

     auto settings = R"( 
        { 
        "geometry" : {"shape" : [128,128, 128] , "domain" : [ [-5,5] , [-5,5] , [-5,5] ] , "coordinates" : "cartesian"} } )"_json;

        
    auto [ba , geom , dm, low_bc, high_bc] = gp::createGeometry(settings["geometry"]);

    amrex::MultiFab phi_new;
    amrex::MultiFab phi_old;

    int nComp=2;
    int nGhosts=2;


    phi_new.define(ba, dm, nComp, nGhosts);
    phi_old.define(ba, dm, nComp, nGhosts);



    gp::wavefunction waveNew(&phi_new,&geom);
    gp::wavefunction waveOld(&phi_old,&geom);
    


    auto it = waveNew.beginAmrexIterator();

    auto region = waveNew[it];

    auto  arr = region.getPhi<gp::array_t>();

    gaussian fillInitial(1.);


    gp::cpuKernelLaunch testKernelLauncher;
    gp::operators::setPhi setPhi;


    testKernelLauncher.apply(waveOld, setPhi,fillInitial);

    testKernelLauncher.apply(waveOld,fillInitial,check_t{1e-6});

    waveOld.fillBoundaries();  


    enum orderDerivative { first = 1};

    auto f = gp::operators::laplacian<orderDerivative::first,gp::DIMENSIONS>{} ;

    testKernelLauncher.apply(waveNew, waveOld, setPhi,f);

    testKernelLauncher.apply(waveNew,fillInitial,lapCheck_t{5e-2});

    gp::trappedGPFunctional<gp::operators::laplacian<1,gp::DIMENSIONS> > func;

    func.evaluate(waveNew,waveOld,0);
    testKernelLauncher.apply(waveNew,fillInitial,gpCheck_t{});

    testKernelLauncher.apply(waveOld, setPhi,fillInitial);
    waveOld.fillBoundaries();  
    gp::normalization::normalize(waveOld,{1});
    testKernelLauncher.apply(waveOld,fillInitial,decltype(fillInitial)::checkNormalization_t{});
    
}   

TEST(wavefunction, harmonicPotential)
{

     auto settings = R"( 
        { 
        "geometry" : {"shape" : [128,128, 128] , "domain" : [ [-5,5] , [-5,5] , [-5,5] ] , "coordinates" : "cartesian"} } )"_json;

        
    auto [ba , geom , dm, low_bc, high_bc] = gp::createGeometry(settings["geometry"]);

    amrex::MultiFab phi_new;
    amrex::MultiFab phi_old;

    int nComp=2;
    int nGhosts=2;

    phi_new.define(ba, dm, nComp, nGhosts);
    phi_old.define(ba, dm, nComp, nGhosts);

    gp::wavefunction waveNew(&phi_new,&geom);
    gp::wavefunction waveOld(&phi_old,&geom);


    gaussian fillInitial(1.);
    gp::operators::setPhi setPhi;
    gp::cpuKernelLaunch testKernelLauncher;


    testKernelLauncher.apply(waveNew, setPhi,fillInitial);
    testKernelLauncher.apply(waveNew,fillInitial,check_t{1e-6});

    waveNew.fillBoundaries();

    enum orderDerivative { first = 1};

     gp::trappedGPFunctional<gp::operators::laplacian<1,gp::DIMENSIONS> > func;


    gp::eulerStepper stepper(&func,true);
    gp::Real time=0;
    gp::Real timeStep= 1e-1 * geom.CellSize()[0] * geom.CellSize()[0]/0.5;
    gp::Real maxTime=1000*timeStep;
    while(time<maxTime)
    {
        std::cout << time << std::endl;

        std::swap(waveOld,waveNew);
        stepper.evolve(waveNew,waveOld,time,timeStep);
        time+=timeStep;
        gp::normalization::normalize(waveNew,{1.});



    }

}




