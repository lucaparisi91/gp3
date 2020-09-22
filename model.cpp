#include "model.h"

void fillMultiFab( MultiFab & fab  , const py::array_t<double> & values )
{
    auto r = values.unchecked<3>();

    for (MFIter mfi(fab); mfi.isValid(); ++mfi)
         {
            const Box& bx = mfi.validbox();

            Array4<Real> const& data = fab[mfi].array();

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
};

model::model(const geometry & geomInfo, int laplacianOrder_, int Ncomp_) : 
laplacianOrder(laplacianOrder_  ) ,Nghost(laplacianOrder_ + 1)  , Ncomp(Ncomp_)
{

    // AMREX_SPACEDIM: number of dimensions
    
    Vector<int> is_periodic(AMREX_SPACEDIM,1);  // periodic in all direction by default

    // make BoxArray and Geometry
    {
      IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
      IntVect dom_hi(AMREX_D_DECL(geomInfo.shape[0]-1, geomInfo.shape[1]-1, geomInfo.shape[2]-1));
      Box domain(dom_lo, dom_hi);

      

      // Initialize the boxarray "ba" from the single box "bx"
      ba.define(domain);
      // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
      //ba.maxSize(max_grid_size);

       // This defines the physical box, [-1,1] in each direction.
        RealBox real_box({AMREX_D_DECL( geomInfo.lower_edges[0],geomInfo.lower_edges[1], geomInfo.lower_edges[2]) },
                         {AMREX_D_DECL( geomInfo.higher_edges[0], geomInfo.higher_edges[1],geomInfo.higher_edges[2] )});

        // This defines a Geometry object
        int coord = 0;
        geom.define(domain,&real_box,coord,is_periodic.data());
    }

    
    
    
    // How Boxes are distrubuted among MPI processes
    dmInit.define(ba);

    phi_real.define(ba, dmInit, Ncomp, Nghost);
    phi_imag.define(ba, dmInit, Ncomp, Nghost);

    phi_real=0.;
    phi_imag=0.;
    


    Vector<BCRec> bc(Ncomp);

    for (int n = 0; n < Ncomp; ++n)
        {
            for (int d=0;d<amrex::SpaceDim ; d++)
            {    
                bc[n].setLo(d, BCType::int_dir);  
                bc[n].setHi(d, BCType::int_dir);
        }
    }

    /*
    Builds the laplacian operator to perform derivatives
    */
    LPInfo info;
    
    //info.setMetricTerm(true);
    info.setMaxCoarseningLevel(0);

    linPoissonReal.define({geom}, {ba}, {dmInit}, info);
    linPoissonImag.define({geom}, {ba}, {dmInit}, info);

    //linPoissonReal.setNComp(nComp);
    
    linPoissonReal.setMaxOrder(laplacianOrder);
    linPoissonImag.setMaxOrder(laplacianOrder);


    std::array<LinOpBCType,AMREX_SPACEDIM> bc_lo;
    std::array<LinOpBCType,AMREX_SPACEDIM> bc_hi;


    for (int d=0;d<amrex::SpaceDim;d++)
        {
        bc_lo[d]=LinOpBCType::Periodic;
        bc_hi[d]=LinOpBCType::Periodic;
        }
    
    linPoissonReal.setDomainBC(bc_lo, bc_hi);
    linPoissonImag.setDomainBC(bc_lo, bc_hi);

    linPoissonReal.setLevelBC(0,nullptr);
    linPoissonImag.setLevelBC(0,nullptr);




}


  void model::fill( const py::array_t<double> & real, const py::array_t<double> & imag  )
  {
      fillMultiFab(phi_real,real);
      fillMultiFab(phi_imag,imag);


  }


