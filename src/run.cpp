#include <algorithm>
#include<array>
#include "src/traits.h"
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




void run(py::array_t<std::complex<Real> > initialCondition , const json_t & settings   )
{
    
    initializer::instance().init();
    auto [ box, geom , dm, low_bc , high_bc] = createGeometry(settings["geometry"]);

    int Ncomp = 1;
    int order = settings["laplacianScheme"]["order"];
    int Nghost = order - 1;


    auto func = initializer::instance().getFunctionalFactory().create(settings["model"]);

    stepper * stepper;


    stepper=new eulerStepper(func, true);

    MultiFab phiOld(box, dm, Ncomp*2, Nghost);
    MultiFab phiNew(box, dm, Ncomp*2, Nghost);
    
    phiOld=0.;
    phiNew=0.;


    std::vector<Real> normalizations;
    
    normalizations = settings["wavefunction"]["normalization"].get<std::vector<Real> >();


    fill(phiNew, initialCondition , geom);


    
    std::string pltfile_init = amrex::Concatenate("out/phi",0,5);

    save(phi_real_old, phi_imag_old, geom , settings, 0 , 0 );
    WriteSingleLevelPlotfile(pltfile_init, phiNew , {"phi"}, geom, 0, 0);

    Real timeStep=settings["stepping"]["timeStep"].get<Real>();
    Real maxTime = settings["stepping"]["maxTime"].get<Real>();
    Real maxTimeBlock = settings["stepping"]["maxTimeBlock"].get<Real>();
    Real time=settings["stepping"]["initialTime"].get<Real>();
    Real lastBlockTime=time;


    std::string outputDir=settings["io"];

    long int iBlock=0;


    while(time <= timeMax)
    {
        lastBlockTime=time;

        while (time < lastBlockTime + maxTimeBlock)
        {   
            // swap old and new solutions
            std::swap(waveNew,waveOld);
            // evolution
            stepper.evolve(waveNew, waveOld,time, timeStep);
            time+=timeStep;
       }
       // output 
       {
            std::cout << "----------------------------------" << std::endl;
            std::cout << "Time: " << time << std::endl; 
            std::cout << "Max Real: "<< currentPhi.max(0) << "; Max Imag: "<< currentPhi.max(1)<<std::endl;
            std::cout << "Min Real: "<< currentPhi.max(0) << "; Min Imag: "<< currentPhi.max(1)<<std::endl;


            std::string pltfile = amrex::Concatenate(outputDir + "/phi_real",i+1,5);

            WriteSingleLevelPlotfile(pltfile, currentPhi, {"phi"}, geom, time, 0);

            //WriteSingleLevelPlotfile(pltfile_imag, phi_imag_old, {"phi"}, geom, time, 0);
       };

    } 

    std::cout << "----------------------------------" << std::endl;
    std::cout << "End at time " << time << std::endl; 
    delete func;
    delete stepper;

    

}
