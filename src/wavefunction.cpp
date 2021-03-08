#include "wavefunction.h"

namespace gp
{

wavefunctionRegion wavefunction::operator[](const amrex::MFIter & mfi) 
{
    return gp::wavefunctionRegion(  mfi.validbox() ,    getPhi().fabPtr(mfi)  , & getGeometry() );
}


constWavefunctionRegion wavefunction::operator[](const amrex::MFIter & mfi)  const
{
    return gp::constWavefunctionRegion(  mfi.validbox() ,    getPhi().fabPtr(mfi)  , & getGeometry() );
}




}