#include <algorithm>
#include<array>
#include "traits.h"
using Real = double;
#include <AMReX_PlotFileUtil.H>
#include "run.h"
#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>
#include "evaluate.h"
#include "timer.h"
#include "stepper.h"
#include "tools.h"
#include <complex>
using namespace amrex;
#include "initializer.h"
#include "plotFile.h"
#include <filesystem>

namespace  fs = std::filesystem ;

inline auto index_F(int i,int j,int xlen,int ylen ) {return i + j*xlen ;}


void fillBox( const Box & bx,  const Array4<Real> & data  , py::array_t<double> values , int ncomps=1)
{
    auto r = values.unchecked<3>();
    const int * lo =bx.loVect();
    const int * hi =bx.hiVect();
    auto size= bx.size();

    int c = 0;
    for (int k=lo[2];k<=hi[2];k++)
        for(int j=lo[1];j<=hi[1];j++)
            for (int i=lo[0];i<=hi[0];i++)
                {
                    data(i,j,k,c)=r(i,j,k);
                }    
}


Real normCylindrical2( const MultiFab & phi_real , const MultiFab & phi_imag,  const Geometry & geom, int component=0)
{
    const Real* dx = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();

    Real norm2=0;
    for ( MFIter mfi( phi_real); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.validbox();
        const int* lo = bx.loVect();
        const int *hi= bx.hiVect();
        Array4< const Real> const & phi_real_box = phi_real[mfi].const_array();
        Array4< const Real> const & phi_imag_box = phi_imag[mfi].const_array();

        for (int j=lo[1];j<=hi[1];j++)
            for (int i=lo[0];i<=hi[0];i++)
            {
                Real r= prob_lo[0] + (i+0.5) * dx[0];

                norm2+=(phi_real_box(i,j,0,component)*phi_real_box(i,j,0,component)
                    + phi_imag_box(i,j,0,component)*phi_imag_box(i,j,0,component))*r;
            }
    }
    amrex::ParallelAllReduce::Sum(norm2,MPI_COMM_WORLD);
    return norm2*2*M_PI*dx[0]*dx[1];
}

void save(MultiFab & realWave, MultiFab & imagWave, Geometry & geom , const json_t & settings, int iBlock , Real time )
    {
        std::string dirname = "out";

        if (settings.contains("name") )
        {
            dirname = settings["name"].get<std::string>();
        }

        if (! fs::exists(dirname)) 
        {
            fs::create_directory(dirname);
        }

        writeSingleLevel(realWave,imagWave,geom,dirname + "/block_" + std::to_string(iBlock),time);

    }

void run(py::array_t<std::complex<Real> > initialCondition , const json_t & settings   )
{
    
    initializer::instance().init();
    auto [ box, geom , dm, low_bc , high_bc] = createGeometry(settings["geometry"]);

    int Ncomp = 1;
    int order = settings["functional"]["laplacian"]["order"];
    int Nghost = order + 1;


    auto func = initializer::instance().getFunctionalFactory().create(settings["functional"]);


    func->define(geom,box,dm, low_bc, high_bc);

    RK4Stepper stepper_phi(func, true,Ncomp,Nghost);

    //eulerStepper stepper_phi(func, true);


    MultiFab phi_real_old(box, dm, Ncomp, Nghost);
    MultiFab phi_imag_old(box, dm, Ncomp, Nghost);

    MultiFab phi_real_new(box, dm, Ncomp, Nghost);
    MultiFab phi_imag_new(box, dm, Ncomp, Nghost);

    phi_imag_new=0.;
    phi_real_new=0.;

    phi_real_old=0.;
    phi_imag_new=0.;


    Real normalization = settings["normalization"].get<Real>();

    fill(phi_real_old, phi_imag_old, initialCondition , geom);
    normalize(phi_real_old,phi_imag_old,geom,normalization);

    
    std::string pltfile_real_init = amrex::Concatenate("out/phi_real",0,5);
    std::string pltfile_imag_init = amrex::Concatenate("out/phi_imag",0,5);


    save(phi_real_old, phi_imag_old, geom , settings, 0 , 0 );


    //WriteSingleLevelPlotfile(pltfile_real_init, phi_real_old, {"phi"}, geom, 0, 0);
    //WriteSingleLevelPlotfile(pltfile_imag_init, phi_imag_old, {"phi"}, geom, 0, 0);


    Real dt=settings["run"]["timeStep"].get<Real>();

    
    
    int nBlocks=settings["run"]["nBlocks"].get<int>() ;
    int stepsPerBlock=settings["run"]["stepsPerBlock"].get<int>() ;
    
    Real time=0;
    for (int i=0;i<nBlocks;i++)
    {
        for (int j=0;j<stepsPerBlock;j++)
        {
            // evolution
            stepper_phi.evolve(phi_real_new,phi_imag_new,phi_real_old,phi_imag_old, time, dt);
            // normalization
            normalize(phi_real_new,phi_imag_new,geom,normalization);
            
            // swap old and new solutions
            std::swap(phi_real_new,phi_real_old);
            std::swap(phi_imag_new,phi_imag_old);

            time+=dt;
       }
       // output 
       {
            std::cout << "----------------------------------" << std::endl;
            std::cout << "Time: " << time << std::endl; 
            std::cout << "Max Real: "<< phi_real_old.max(0) << "; Max Imag: "<< phi_imag_old.max(0)<<std::endl;
            std::cout << "Min Real: "<< phi_real_old.min(0) << "; Min Imag: "<< phi_imag_old.min(0)<<std::endl;

            std::string pltfile_real = amrex::Concatenate("out/phi_real",i+1,5);
            std::string pltfile_imag = amrex::Concatenate("out/phi_imag",i+1,5);

            save(phi_real_old, phi_imag_old, geom , settings, i , time );

            //WriteSingleLevelPlotfile(pltfile_real, phi_real_old, {"phi"}, geom, time, 0);
            //WriteSingleLevelPlotfile(pltfile_imag, phi_imag_old, {"phi"}, geom, time, 0);
       }
    } 

    



    std::cout << "----------------------------------" << std::endl;
    std::cout << "End at time " << time << std::endl; 
    delete func;

}
