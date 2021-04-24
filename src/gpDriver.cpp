#include "traits.h"
#include "wavefunction.h"
#include "gpDriver.h"
#include "normalization.h"
#include "functional.h"
#include "geometry.h"
#include "stepper.h"
#include <string>
#include "timers.h"
#include "functionalFactory.h"

namespace gp
{
namespace gpDriver
{


using json_t = gp::json_t;
using Real = gp::Real;


void run(const json_t & settings)
{
    auto [ geom , low_bc, high_bc] = 

    gp::createGeometry(settings["geometry"]);

     

    auto  phi_new = gp::createMultiFab(settings);
    auto phi_old = gp::createMultiFab(settings);
    
    gp::Real N = settings["normalization"];


    // read initial condition
    json_t jWave;
    std::ifstream f;
    
    f.open(settings["initialCondition"].get<std::string>() + std::string("/wave.json"));
    f >> jWave ;
    f.close();



    auto  phi0 = gp::createMultiFab(jWave);


    phi_new.ParallelCopy(phi0);
    phi_old.ParallelCopy(phi0);




    gp::wavefunction waveNew(&phi_new,&geom);
    gp::wavefunction waveOld(&phi_old,&geom);


    waveNew.fillBoundaries();


    // create functional

    gp::normalization::normalize(waveNew,{N});

    
    gp::functionalFactory fac;
    fac.registerFunctional<gp::trappedGPFunctional<gp::operators::laplacian<1,gp::DIMENSIONS> > >("trappedGP");

    auto func = fac.create(settings["functional"]);



    bool imaginaryTime=settings["imaginaryTime"].get<bool>();

    gp::eulerStepper stepper(&(*func),imaginaryTime);
    gp::Real time=settings["initialTime"].get<gp::Real>();
    gp::Real timeStep=settings["timeStep"].get<gp::Real>();


    gp::Real maxTime = settings["maxTime"].get<gp::Real>();


    size_t iteration=0;
    size_t stepsPerBlock= settings["stepsPerBlock"].get<size_t>();

    // save the initial condition
    waveNew.save("phi" + std::to_string(iteration));

    START_TIMER("timeEvolution");

    while(time<maxTime)
    {
        std::swap(waveOld,waveNew);
        START_TIMER("stepping");
        stepper.evolve(waveNew,waveOld,time,timeStep);
        STOP_TIMER("stepping");
        time+=timeStep;
        START_TIMER("normalize");

        gp::normalization::normalize(waveNew,{1.});
        STOP_TIMER("normalize");
        
        iteration+=1;

        START_TIMER("out");
           if ( (iteration ) % stepsPerBlock == 0)
        {

        START_TIMER("max");
        gp::Real max = waveNew.getPhi().max(0);
        STOP_TIMER("max");


        START_TIMER("save");
        waveNew.save("phi" + std::to_string(iteration));
        STOP_TIMER("save");
            if (amrex::ParallelDescriptor::IOProcessor())
        {
            std::cout << "Time: " << time << std::endl;
            std::cout << "Max: " << max << std::endl;
        }

       }
       STOP_TIMER("out");    
     
    }
    STOP_TIMER("timeEvolution");

    
    auto report = timers::getInstance().report();
    if (amrex::ParallelDescriptor::IOProcessor())
    {
        std::cout << report << std::endl ;
    }
   
}



}
}