#include "functional.h"

namespace gp {


 template<class laplacianOperator_t>
void trappedGPFunctional<laplacianOperator_t>::evaluate(
	wavefunction & waveNew,const  wavefunction & waveOld, Real time )
{

    const auto & geom = waveNew.getGeometry();

    for ( auto mfi = waveNew.beginAmrexIterator() ; mfi.isValid(); ++mfi ) 
        { 
            auto  waveNewRegion = waveNew[mfi];
            auto  waveOldRegion = waveOld[mfi];

            const auto & phiNew = waveNewRegion.getPhi<gp::array_t>();
            const auto & phiOld = waveOldRegion.getPhi<gp::array_t>();

        if constexpr ( DIMENSIONS == 3)
            for (int k=waveNewRegion.minIndex(2);k<=waveNewRegion.maxIndex(2);k++)
                for (int j=waveNewRegion.minIndex(1);j<=waveNewRegion.maxIndex(1);j++) 
                    for (int i=waveNewRegion.minIndex(0);i<=waveNewRegion.maxIndex(0);i++)
                    {
                        
                    gp::Real x = geom.ProbLo()[0] + geom.CellSize()[0]*i  ;
                    gp::Real y = geom.ProbLo()[1] + geom.CellSize()[1]*j  ;
                    gp::Real z = geom.ProbLo()[2] + geom.CellSize()[2]*k  ;


                    for(int c=0;c<2*nComp;c++)
                    {
                        auto tmp = -0.5*_lap(i,j,k,c,phiOld,geom) + prefactor*(x*x + y*y + z*z)
                        ;

                        phiNew(i,j,k,c)=tmp;
                    }


                    }

    }
};





template class trappedGPFunctional<operators::laplacian<1,DIMENSIONS> >;


}