#include "gtest/gtest.h"
#include "geometry.h"
#include "wavefunction.h"
#include "kernelLauncher.h"


struct check_t{
    gp::Real tol = 1e-5;
};

struct lapCheck_t {gp::Real tol = 1e-4; };


struct gaussian : public gp::operators::operatorBase
{

    gaussian(const gp::Real alpha_) : alpha(alpha_)
    {

    };

    gp::Real operator()(int i, int j , int k,gp::wavefunctionRegion & wave, int c)
    {
        auto & phi = wave.getPhiNew<gp::array_t>();

        const auto & geom = wave.getGeometry();

        gp::Real x = geom.ProbLo()[0] + geom.CellSize()[0]*i  ;
        gp::Real y = geom.ProbLo()[1] + geom.CellSize()[1]*j  ;
        gp::Real z = geom.ProbLo()[2] + geom.CellSize()[2]*k  ;


        return exp(-alpha*(x*x + y*y + z*z)  );
           
    }

    void operator()(int i, int j , int k, gp::wavefunctionRegion & wave, const check_t & check)
    {
        auto & phi = wave.getPhiNew<gp::array_t>();

        const auto & geom = wave.getGeometry();

        gp::Real x = geom.ProbLo()[0] + geom.CellSize()[0]*i  ;
        gp::Real y = geom.ProbLo()[1] + geom.CellSize()[1]*j  ;
        gp::Real z = geom.ProbLo()[2] + geom.CellSize()[2]*k  ;

        ASSERT_NEAR( phi(i,j,k)  , exp(-alpha*(x*x + y*y + z*z))  , check.tol );

    }

    void operator()(int i, int j , int k, gp::wavefunctionRegion & wave, const lapCheck_t & check)
    {
        auto & phi = wave.getPhiNew<gp::array_t>();

        const auto & geom = wave.getGeometry();

        gp::Real x = geom.ProbLo()[0] + geom.CellSize()[0]*i  ;
        gp::Real y = geom.ProbLo()[1] + geom.CellSize()[1]*j  ;
        gp::Real z = geom.ProbLo()[2] + geom.CellSize()[2]*k  ;
        auto r2 = x*x + y*y + z*z;

        ASSERT_NEAR( phi(i,j,k)  , exp(-alpha*r2) *(-2*alpha*(3-2*alpha*r2) )  , check.tol );

    }

    private:
    const gp::Real alpha;

};

TEST(initialization, createGeometry)
{


    auto settings = R"( 
        { 
        "geometry" : {"shape" : [200,200, 200] , "domain" : [ [-5,5] , [-5,5] , [-5,5] ] , "coordinates" : "cartesian"} } )"_json;
    
    auto [box , geom , dm, low_bc, high_bc] = gp::createGeometry(settings["geometry"]);

    int Ncomp = 2;
    int Nghost = 2;

    gp::multiFab phi_new(box, dm, Ncomp, Nghost);
    gp::multiFab phi_old(box, dm, Ncomp, Nghost);


    gp::wavefunction wave(&phi_new,&phi_old,&geom);

    auto it = wave.beginAmrexIterator();

    auto region = wave[it];

    auto  arr = region.getPhiNew<gp::array_t>();

    gaussian fillInitial(1.);


    gp::cpuKernelLaunch testKernelLauncher;
    gp::operators::setPhiNew setPhiNew;


    testKernelLauncher.apply(wave,setPhiNew,fillInitial);

    testKernelLauncher.apply(wave,fillInitial,check_t{1e-6});

    wave.fillBoundaries();
    wave.swapOldAndNew();   


    enum orderDerivative { first = 1};

    auto f = gp::operators::laplacian<orderDerivative::first,gp::DIMENSIONS>{} ; // defines the function of the PDE dphidt = f(phi)


    testKernelLauncher.apply(wave, setPhiNew,f) ;

    testKernelLauncher.apply(wave,fillInitial,lapCheck_t{5e-2});
    
    
}