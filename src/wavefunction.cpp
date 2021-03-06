#include "wavefunction.h"

namespace gp
{

wavefunctionRegion wavefunction::operator[](const amrex::MFIter & mfi)
{
    return gp::wavefunctionRegion(  mfi.validbox() ,    getPhiNew().fabPtr(mfi) , getPhiOld().fabPtr(mfi) , & getGeometry() );
}


}