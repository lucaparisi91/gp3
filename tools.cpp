#include "tools.h"

initializer *initializer::s_instance = 0;

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