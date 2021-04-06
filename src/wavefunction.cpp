#include "wavefunction.h"
#include <AMReX_VisMF.H>
#include <filesystem>

namespace fs = std::filesystem;

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


void wavefunction::save(const std::string & folder)
{

     auto amrexFileName = folder + std::string("/") + _name ;
    
    
    // save multifab
    if (amrex::ParallelDescriptor::IOProcessor() )
    {
         if (!fs::exists(folder)) 
        {
            fs::create_directory(folder);
        } 

    } 

    amrex::ParallelDescriptor::Barrier();

    amrex::VisMF::Write(*_phi, amrexFileName);

    if (amrex::ParallelDescriptor::IOProcessor() )
     {
    auto jGeo= toJson( getGeometry());


    json_t j;

    const auto & ba = (*_phi).boxArray();
    auto jsonFilename = folder + std::string("/") + std::string("wave.json") ;

    j["boxes"]=toJson(ba);
    j["geometry"]=jGeo;
    j["components"]=_phi->nComp()/2;

    std::vector<int> nGhosts;

    for(int d=0;d<DIMENSIONS;d++)
    {
        nGhosts.push_back(_phi->nGrow(d) );
    }
    
    j["nGhosts"]=nGhosts;
    j["name"]=_name;
    j["folder"]=folder;

    std::ofstream f;
    f.open(jsonFilename);

    f << j ;

    f.close();

     }



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

    //std::cout << "n boxes: " << ba.size() << std::endl;


    amrex::DistributionMapping dm(ba);
    
    phi.define(ba, dm, nComp, nGhostsAmrex); 
    phi=0;
    

    if ( 
        (settings.find("folder") !=settings.end() ) and 
        (settings.find("name") != settings.end() )
        )
        {
            auto filename = settings["folder"].get<std::string>() + std::string("/") + settings["name"].get<std::string>() ;
            std::cout << "filename: " << filename << std::endl;
            if (fs::exists(filename + "_H")) {
            amrex::VisMF::Read(phi, filename);
            }
        }


    



    return phi;
}





}