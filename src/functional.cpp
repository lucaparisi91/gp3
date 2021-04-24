#include "functional.h"
#include <omp.h>

namespace gp {

 template<class laplacianOperator_t>
void trappedGPFunctional<laplacianOperator_t>::evaluate(
	wavefunction & waveNew,const  wavefunction & waveOld, Real time )
{
    const auto & geom = waveNew.getGeometry();
    
    _lap.setGeometry(geom);


    auto & mfNew = waveNew.getPhi();
    const auto & mfOld = waveOld.getPhi();
{
#pragma omp parallel
    for (amrex::MFIter mfi(mfNew,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<Real> const& phiNew = mfNew.array(mfi); 
            amrex::Array4<const Real> const& phiOld = mfOld.const_array(mfi); 

            
            int i0=bx.smallEnd(0); int i1=bx.bigEnd(0);

            const int C= waveNew.nComp();

            std::array<gp::Real,3> probLo={ AMREX_D_DECL(geom.ProbLo()[0],geom.ProbLo()[1],geom.ProbLo()[2])};
            std::array<gp::Real,3> cellSize={ AMREX_D_DECL(geom.CellSize()[0],geom.CellSize()[1],geom.CellSize()[2])};

        if constexpr ( DIMENSIONS == 3)
        {
#pragma omp for
            for (int k=bx.smallEnd(2);k<=bx.bigEnd(2);k++)
            {
                for (int j=bx.smallEnd(1);j<=bx.bigEnd(1);j++) 
                    #pragma omp simd
                    for (int i=i0;i<=i1;i++)
                    {
                
                    gp::Real x = probLo[0] + cellSize[0]*(i+0.5)  ;
                    gp::Real y = probLo[1] + cellSize[1]*(j+0.5)  ;
                    gp::Real z = probLo[2] + cellSize[2]*(k + 0.5)  ;

                    for(int c=0;c<C;c++)
                    {
                        phiNew(i,j,k,c) = -0.5*_lap(i,j,k,c,phiOld) + prefactor*(x*x + y*y + z*z)*phiOld(i,j,k,c);
                        
                    }
                    
                    }

            }

        }
    }
}
};

 template<class laplacianOperator_t>
void dropletFunctional<laplacianOperator_t>::evaluate(
	wavefunction & waveNew,const  wavefunction & waveOld, Real time )
{
    const auto & geom = waveNew.getGeometry();

    
    _lap.setGeometry(geom);

    auto & mfNew = waveNew.getPhi();
    const auto & mfOld = waveOld.getPhi();
{
#pragma omp parallel
    for (amrex::MFIter mfi(mfNew,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<Real> const& phiNew = mfNew.array(mfi); 
            amrex::Array4<const Real> const& phiOld = mfOld.const_array(mfi); 

            int i0=bx.smallEnd(0); int i1=bx.bigEnd(0);

            const int C= waveNew.nComp();

            assert(C==2);


            std::array<gp::Real,3> probLo={ AMREX_D_DECL(geom.ProbLo()[0],geom.ProbLo()[1],geom.ProbLo()[2])};
            std::array<gp::Real,3> cellSize={ AMREX_D_DECL(geom.CellSize()[0],geom.CellSize()[1],geom.CellSize()[2])};


        if constexpr ( DIMENSIONS == 3)
        {
#pragma omp for
            for (int k=bx.smallEnd(2);k<=bx.bigEnd(2);k++)
            {
                for (int j=bx.smallEnd(1);j<=bx.bigEnd(1);j++) 
                    #pragma omp simd
                    for (int i=i0;i<=i1;i++)
                    {

                        auto density=phiOld(i,j,k,0)*phiOld(i,j,k,0) +
                        phiOld(i,j,k,1)*phiOld(i,j,k,1);

                        auto tmp=  -3*density
                        + 5/2. * std::pow(density,3/2.);

                        phiNew(i,j,k,0)=-0.5*_lap(i,j,k,0,phiOld) + tmp *phiOld(i,j,k,0);
                        phiNew(i,j,k,1)=-0.5*_lap(i,j,k,1,phiOld) + tmp *phiOld(i,j,k,1);

                    
                    }

            }

        }
    }
}
};

void setGaussian(
	wavefunction & waveNew, Real alpha, int c)
{

    const auto & geom = waveNew.getGeometry();

    for ( auto mfi = waveNew.beginAmrexIterator() ; mfi.isValid(); ++mfi ) 
        { 
            auto  waveNewRegion = waveNew[mfi];
           
            const auto & phiNew = waveNewRegion.getPhi<gp::array_t>();
           
        if constexpr ( DIMENSIONS == 3)
            for (int k=waveNewRegion.minIndex(2);k<=waveNewRegion.maxIndex(2);k++)
                for (int j=waveNewRegion.minIndex(1);j<=waveNewRegion.maxIndex(1);j++) 
                    for (int i=waveNewRegion.minIndex(0);i<=waveNewRegion.maxIndex(0);i++)
                    {
                        
                    gp::Real x = geom.ProbLo()[0] + geom.CellSize()[0]*(i+0.5)  ;
                    gp::Real y = geom.ProbLo()[1] + geom.CellSize()[1]*(j+0.5)  ;
                    gp::Real z = geom.ProbLo()[2] + geom.CellSize()[2]*(k+0.5)  ;


                    phiNew(i,j,k,c)=std::exp(-alpha*(x*x + y*y + z*z));
                    }


        }

}

template class trappedGPFunctional<operators::laplacian<1,DIMENSIONS> >;


}