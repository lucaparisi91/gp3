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

using namespace amrex;

inline auto index_F(int i,int j,int xlen,int ylen ) {return i + j*xlen ;} 

void fillBox( const Box & bx,  const Array4<Real> & data  , py::array_t<double> values )
{
    auto r = values.unchecked<2>();  
    const int * lo =bx.loVect();
    const int * hi =bx.hiVect();
    auto size= bx.size();

    for(int j=lo[1];j<=hi[1];j++)
        for (int i=lo[0];i<=hi[0];i++)
            {
                data(i,j,0)=r(i,j);
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

void normalize(  MultiFab & phi_real ,  MultiFab & phi_imag,  const Geometry & geom)
{
    /* Normalizes each component indipendently  */
    for (int i=0;i<phi_real.nComp();i++)
        {
            Real norm2=normCylindrical2(phi_real,phi_imag,geom,i);
            phi_real.mult(1./std::sqrt(norm2),i,1);
            phi_imag.mult(1./std::sqrt(norm2),i,1);
    }
}


void run(py::array_t<double> initialCondition_real,py::array_t<double> initialCondition_imag,const geometry & geomInfo)
{	
    amrex::Initialize(MPI_COMM_WORLD);

	  // What time is it now?  We'll use this to compute total run time.
    Real strt_time = amrex::second();
    
    // AMREX_SPACEDIM: number of dimensions
    int n_cell, max_grid_size, maxSteps, plot_int;
    Vector<int> is_periodic(AMREX_SPACEDIM,1);  // periodic in all direction by default

    is_periodic[0]=0;
    // make BoxArray and Geometry
    BoxArray ba;
    Geometry geom;
    {
      IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
      IntVect dom_hi(AMREX_D_DECL(geomInfo.shape[0]-1, geomInfo.shape[1]-1, 0));
      Box domain(dom_lo, dom_hi);

        // Initialize the boxarray "ba" from the single box "bx"
        ba.define(domain);
        // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
        //ba.maxSize(max_grid_size);

       // This defines the physical box, [-1,1] in each direction.
        RealBox real_box({AMREX_D_DECL( geomInfo.lower_edges[0],geomInfo.lower_edges[1], 0.) },
                         {AMREX_D_DECL( geomInfo.higher_edges[0], geomInfo.higher_edges[1], 0.)});

        // This defines a Geometry object
        geom.define(domain,&real_box,1,is_periodic.data());
    }

    // Nghost = number of ghost cells for each array 
    int Nghost = 3;
    
    // Ncomp = number of components for each array
    int Ncomp  = 1;
    
    // How Boxes are distrubuted among MPI processes
    DistributionMapping dmInit(ba);
    
    // we allocate two phi multifabs; one will store the old state, the other the new.

    MultiFab phi_real_old(ba, dmInit, Ncomp, Nghost);
    MultiFab phi_imag_old(ba, dmInit, Ncomp, Nghost);

    MultiFab phi_real_new(ba, dmInit, Ncomp, Nghost);
    MultiFab phi_imag_new(ba, dmInit, Ncomp, Nghost);

    phi_imag_new=0.;
    phi_real_new=0.;

    phi_real_old=0.;
    phi_imag_new=0.;    

    for (MFIter mfi(phi_real_new); mfi.isValid(); ++mfi)
         {
            const Box& bx = mfi.validbox();

            Array4<Real> const& data_real = phi_real_old[mfi].array();
            Array4<Real> const& data_imag = phi_imag_old[mfi].array();

            fillBox(bx,data_real,initialCondition_real);
            fillBox(bx,data_imag,initialCondition_imag);

        }

    std::string pltfile_real_init = amrex::Concatenate("out/phi_real",0,5);
    std::string pltfile_imag_init = amrex::Concatenate("out/phi_imag",0,5);

    Vector<BCRec> bc(Ncomp);
    for (int n = 0; n < Ncomp; ++n)
        {
            bc[n].setLo(0, BCType::foextrap); // first-order extrapolation     
            bc[n].setLo(1, BCType::int_dir);  // external Dirichlet

            bc[n].setHi(0, BCType::foextrap); // first-order extrapolation     
            bc[n].setHi(1, BCType::int_dir);  // external Dirichlet
        }
    /*
    Builds the laplacian operator to perform derivatives
    */
    LPInfo info;


    info.setMetricTerm(true);
    info.setMaxCoarseningLevel(0);

    MLPoisson linPoissonReal({geom}, {ba}, {dmInit}, info);
    MLPoisson linPoissonImag({geom}, {ba}, {dmInit}, info);

    //linPoissonReal.setNComp(nComp);

    // order of stencil
    int linop_maxorder = 2;
    linPoissonReal.setMaxOrder(linop_maxorder);
    linPoissonImag.setMaxOrder(linop_maxorder);

    std::array<LinOpBCType,AMREX_SPACEDIM> bc_lo;
    std::array<LinOpBCType,AMREX_SPACEDIM> bc_hi;

    
    bc_lo[1]=LinOpBCType::Periodic;
    bc_hi[1]=LinOpBCType::Periodic;

    bc_lo[0]=LinOpBCType::Neumann;
    bc_hi[0]=LinOpBCType::Neumann;

    linPoissonReal.setDomainBC(bc_lo, bc_hi);
    linPoissonImag.setDomainBC(bc_lo, bc_hi);

    RK4Stepper stepper_phi(false,Ncomp,Nghost);
    //eulerStepper stepper_phi(true);
    const Real* dx = geom.CellSize();


    normalize(phi_real_old,phi_imag_old,geom);

    WriteSingleLevelPlotfile(pltfile_real_init, phi_real_old, {"phi"}, geom, 0, 0);
    WriteSingleLevelPlotfile(pltfile_imag_init, phi_imag_old, {"phi"}, geom, 0, 0);

    Real dt=1e-5;
    linPoissonReal.setLevelBC(0,nullptr);
    linPoissonImag.setLevelBC(0,nullptr);

    int nBlocks=10000;
    int stepsPerBlock=1000;

    Real time=0;
    for (int i=0;i<nBlocks;i++)
    {
        for (int j=0;j<stepsPerBlock;j++)
        {
            // evolution
            stepper_phi.evolve(phi_real_new,phi_imag_new,phi_real_old,phi_imag_old,linPoissonReal,linPoissonImag,geom, bc, time, dt);
            // normalization
            normalize(phi_real_new,phi_imag_new,geom);
            
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

            WriteSingleLevelPlotfile(pltfile_real, phi_real_old, {"phi"}, geom, time, 0);
            WriteSingleLevelPlotfile(pltfile_imag, phi_imag_old, {"phi"}, geom, time, 0);

       }
    }   

    std::cout << "----------------------------------" << std::endl;
    std::cout << "End at time " << time << std::endl; 

}