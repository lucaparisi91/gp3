#include "normalization.h"

namespace gp
{

    namespace normalization
    {
        std::vector<Real> normRealWavefunctionCartesian(const wavefunction & wave)
        {
             std::vector<Real> currentNorms(wave.nComp() ,0);
             return currentNorms;
        }

        std::vector<Real> normComplexWavefunctionCartesian(const wavefunction & wave)
        {
            assert(wave.isReal() == false);
            assert(wave.nComp() % 2 == 0);


            const auto & phi = wave.getPhi();
            const auto & geom = wave.getGeometry();

            std::vector<Real> currentNorms(wave.nComp()/2);

            Real dV=1;
            for(int d=0;d<DIMENSIONS;d++)
            {
                dV*=geom.CellSize()[d];
            }

            for (int c=0;c<wave.nComp()/2;c+=1)
            {
                Real tmp,tmp2;
                tmp=phi.norm2(2*c)*std::sqrt(dV);
                tmp2=phi.norm2(2*c+1)*std::sqrt(dV);
                
                currentNorms[c]=std::sqrt(tmp*tmp + tmp2*tmp2);

                //amrex::ParallelAllReduce::Sum(currentNorms[c],MPI_COMM_WORLD);


            }

            return currentNorms;
        }


        std::vector<Real> norm(const wavefunction & wave)
        {
            if ( not wave.isReal() )
            {
                return normComplexWavefunctionCartesian(wave);
            }
            else
            {
                return normRealWavefunctionCartesian(wave);
            }


        }


        void normalize(wavefunction & wave,const std::vector<Real> & norms)
        {
            auto currentNorms=norm(wave);

            auto & phi = wave.getPhi();

            if (not wave.isReal())
            {
                for (int c=0;c<wave.nComp()/2;c+=1)
                {
                    Real C = norms[c] /currentNorms[c];

                    phi.mult(C, c, 2, phi.nGrow() );

                }
            }

        }



    }

};