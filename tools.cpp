#include "tools.h"
#include "gpExceptions.h"

Real norm( const MultiFab & phi_real , const MultiFab & phi_imag,  const Geometry & geom, int component)
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

#if AMREX_SPACEDIM == 2
        for (int j=lo[1];j<=hi[1];j++)
            for (int i=lo[0];i<=hi[0];i++)
            {
             
                norm2+=(phi_real_box(i,j,0,component)*phi_real_box(i,j,0,component)
                    + phi_imag_box(i,j,0,component)*phi_imag_box(i,j,0,component));
            }
    
#endif


#if AMREX_SPACEDIM == 3
        for (int k=lo[2];k<=hi[2];k++)
            for (int j=lo[1];j<=hi[1];j++)
                for (int i=lo[0];i<=hi[0];i++)
                    {

                        norm2+=(phi_real_box(i,j,k,component)*phi_real_box(i,j,k,component)
                        + phi_imag_box(i,j,k,component)*phi_imag_box(i,j,k,component));
                    }
#endif      
    }

    amrex::ParallelAllReduce::Sum(norm2,MPI_COMM_WORLD);


    Real dV=1;

    for (int d=0;d<amrex::SpaceDim;d++)
    {
        dV*=dx[d];
    }

    return norm2*dV;
}





std::tuple< BoxArray , Geometry , DistributionMapping   >  
createGeometry( const json_t & settings)
{
    BoxArray ba;
    Geometry geom;
    Vector<int> is_periodic(AMREX_SPACEDIM,1);

    std::array<size_t,AMREX_SPACEDIM> shape;
    std::array<Real,AMREX_SPACEDIM> lower_edges;
    std::array<Real,AMREX_SPACEDIM> higher_edges;
    
    std::string coordinates = settings["coordinates"].get<std::string>();


    for (int i=0;i<AMREX_SPACEDIM;i++)
    {
        lower_edges[i]=settings["domain"][i][0].get<Real>();
        higher_edges[i]=settings["domain"][i][1].get<Real>();
        shape[i] = settings["shape"][i].get<int>();

    }

    IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
    IntVect dom_hi(AMREX_D_DECL( shape[0]-1, shape[1]-1, shape[2]-1));
    Box domain(dom_lo, dom_hi);
    ba.define(domain);
   //ba.maxSize(max_grid_size);

    RealBox real_box({AMREX_D_DECL( lower_edges[0],lower_edges[1],lower_edges[2]) },
                         {AMREX_D_DECL( higher_edges[0], higher_edges[1], higher_edges[2] )});
                        
    int coord = 0;

    if (coordinates == "cartesian")
    {
        coord=0;
    }
    else if (coordinates == "spherical")
    {
        coord=2;

        is_periodic[0]=0;


    }
    else
    {
        throw missingImplementation("Unkown coordinates :" + coordinates );
    }

    geom.define(domain,&real_box,coord,is_periodic.data());

    DistributionMapping dm(ba);

    return {ba, geom, dm} ;

}

void fill(MultiFab & realState, MultiFab & imagState, py::array_t<std::complex<Real> > initialCondition , Geometry & geom)
{
    auto psi = initialCondition.unchecked<AMREX_SPACEDIM>();

    LOOP(realState,geom)

#if AMREX_SPACEDIM == 3
       data(i,j,k)=std::real( psi(i,j,k) );
#endif

#if AMREX_SPACEDIM == 1
       data(i,j,k)=std::real( psi(i) );
#endif



    ENDLOOP

    LOOP(imagState,geom)

#if AMREX_SPACEDIM == 3
       data(i,j,k)=std::imag( psi(i,j,k) );
#endif

#if AMREX_SPACEDIM == 1
       data(i,j,k)=std::imag( psi(i) );
#endif

    ENDLOOP



}
