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


void run(py::array_t<double> initialCondition,const geometry & geomInfo)
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

            Array4<Real> const& data = phi_real_old[mfi].array();
            fillBox(bx,data,initialCondition);
        }

    std::string pltfile = amrex::Concatenate("out/plt",0,5);

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

    eulerStepper stepper_phi;

    WriteSingleLevelPlotfile(pltfile, phi_real_old, {"phi"}, geom, 0, 0);
    Real dt=0.01;
    

    stepper_phi.evolve(phi_real_new,phi_imag_new,phi_real_old,phi_imag_old,linPoissonReal,linPoissonImag,geom, bc, 0., dt);
    

    pltfile = amrex::Concatenate("out/plt",1,5);
    WriteSingleLevelPlotfile(pltfile, phi_real_new, {"phi"}, geom, 0, 0);

}