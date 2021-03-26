#include "wavefunction.h"
#include <AMReX_VisMF.H>




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




amrex::MultiFab  createMultiFab(const json_t & settings)
{
    auto [geom , low_bc, high_bc] = gp::createGeometry(settings["geometry"]);

      // create multifab
     amrex::MultiFab phi;
    int nComp = 2;


    auto nGhosts= settings["nGhosts"].get<std::vector<int> >();

    amrex::IntVect nGhostsAmrex(AMREX_D_DECL(nGhosts[0],nGhosts[1],nGhosts[2]));


    if (settings.find("components") != settings.end())
        {
            nComp=settings["components"].get<int>()*2;

        }
 
    auto ba = createGrids(settings);

    amrex::DistributionMapping dm(ba);


    phi.define(ba, dm, nComp, nGhostsAmrex); 


    if ( 
        (settings.find("folder") !=settings.end() ) and 
        (settings.find("name") != settings.end() )
        )
        {
            auto filename = settings["folder"].get<std::string>() + std::string("/") + settings["name"].get<std::string>() ;


            amrex::VisMF::Read(phi, filename);
        }
    {

    }


    return phi;
}





}