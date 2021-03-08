#include "operators.h"

namespace gp
{

class cpuKernelLaunch
{

public:
    cpuKernelLaunch(){}


    template<class op_t,class ... Args>
    auto apply(wavefunction &  wave, op_t && currentOp ,Args&&...  args)
    {
        applyNoReturn( wave, std::forward<op_t> (currentOp) ,std::forward<Args>(args)...  );
    }
    template<class op_t,class ... Args>
    auto apply(wavefunctionRegion &  wave, op_t && currentOp ,Args&&...  args)
    {
        applyNoReturn( wave, std::forward<op_t> (currentOp) ,std::forward<Args>(args)...  );
    }



    template<class op_t,class ... Args>
    auto applyNoReturn(wavefunction &  wave, op_t && currentOp ,Args&&...  args)
    {
        for ( auto mfi = wave.beginAmrexIterator() ; mfi.isValid(); ++mfi ) 
        { 
            auto  currentWaveRegion = wave[mfi];
            applyNoReturn(currentWaveRegion, std::forward<op_t>(currentOp) , std::forward<Args>(args)...);
        }
        
    };

    template<class op_t,class ... Args>
    auto applyNoReturn(wavefunction &  wave, wavefunction & wave2 ,op_t && currentOp ,Args&&...  args)
    {
        for ( auto mfi = wave.beginAmrexIterator() ; mfi.isValid(); ++mfi ) 
        { 
            auto  currentWaveRegion = wave[mfi];
            auto  currentWaveRegion2 = wave2[mfi];
            
            applyNoReturn(currentWaveRegion,currentWaveRegion2, std::forward<op_t>(currentOp) , std::forward<Args>(args)...);
        }
        
    };



template<class op_t,class ... Args>
    auto applyNoReturn(wavefunctionRegion &  wave, op_t && currentOp ,Args&&...  args)
    {

    if constexpr ( DIMENSIONS == 3)
        for (int k=wave.minIndex(2);k<=wave.maxIndex(2);k++)
	 		for (int j=wave.minIndex(1);j<=wave.maxIndex(1);j++) 
	     		for (int i=wave.minIndex(0);i<=wave.maxIndex(0);i++)
                 {
                    for(int c=0;c<wave.nComp();c++)
                    {
                        currentOp(i,j,k, c, wave,std::forward<Args>(args)... );
                    }

                 }

    }


template<class op_t,class ... Args>
    auto applyNoReturn(wavefunctionRegion &  wave,wavefunctionRegion &  wave2, op_t && currentOp ,Args&&...  args)
    {

    if constexpr ( DIMENSIONS == 3)
        for (int k=wave.minIndex(2);k<=wave.maxIndex(2);k++)
	 		for (int j=wave.minIndex(1);j<=wave.maxIndex(1);j++) 
	     		for (int i=wave.minIndex(0);i<=wave.maxIndex(0);i++)
                 {
                    for(int c=0;c<wave.nComp();c++)
                    {
                    currentOp(i,j,k, c, wave,wave2,std::forward<Args>(args)... );
                    }
                 }

    }
};


}