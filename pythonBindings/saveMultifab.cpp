#include "saveMultifab.h"
#include "traits.h"
#include <utility>
#include <AMReX_BoxArray.H>
#include "geometry.h"
#include "wavefunction.h"
#include <AMReX_VisMF.H>
#include "initializer.h"


auto  get(py::list list, int i)
{
    return *(list.begin() + i); 
}

auto createBoxes(const json_t & j)
{
    // load boxes from json file
    amrex::Vector<amrex::Box> boxes;

    for (auto currentBox : j )
    {

        amrex::IntVect lowerEdge;
        amrex::IntVect upperEdge;

        for (int d=0;d<gp::DIMENSIONS;d++)
        {
            auto limitsDimension=currentBox["range"][d].get<std::pair<int,int> >();
            lowerEdge[d]=limitsDimension.first;
            upperEdge[d]=limitsDimension.second;

        }

        boxes.push_back(amrex::Box(lowerEdge,upperEdge));


    }

    amrex::BoxList bl(std::move(boxes) );
    amrex::BoxArray ba(bl);

    return ba;
}

auto createMultifab(const json_t & settings)
{
    auto ba = createBoxes(settings["boxes"]);

    // create geometry from json file

     auto [baDomain , geom , dm, low_bc, high_bc] = gp::createGeometry(settings["geometry"]);


     // create multifab
     amrex::MultiFab phi;


    auto nGhosts= settings["nGhosts"].get<std::vector<int> >();

    amrex::IntVect nGhostsAmrex(AMREX_D_DECL(nGhosts[0],nGhosts[1],nGhosts[2]));


    int nComp=settings["components"].get<int>()*2;

    phi.define(ba, dm, nComp, nGhostsAmrex); 
    return phi;
}


void saveMultifab( py::list initialConditions , const json_t & settings   )
{
    initializer::getInstance().init();

    auto [baDomain , geom , dm, low_bc, high_bc] = gp::createGeometry(settings["geometry"]);

    auto phi = createMultifab(settings);

    gp::wavefunction wave(&phi,&geom);

    int i=0;
    for ( auto mfi = wave.beginAmrexIterator() ; mfi.isValid(); ++mfi )
    {
        // lop box region
        auto region = wave[mfi];
        auto & phi_arr = region.getPhi<gp::array_t>();
        //std::cout << region.maxIndex(0) << std::endl;
        
        // load python array of data
        auto arr=get(initialConditions,i);

        using pythonFortranArray_t = py::array_t<std::complex<double>,py::array::f_style | py::array::forcecast > ;


        auto boxData  = py::cast<pythonFortranArray_t>(arr);


        py::buffer_info buf1 = boxData.request();

        if (buf1.ndim != ( gp::DIMENSIONS + 1 ) )
        {
            throw std::runtime_error("Number of dimensions must be " + std::to_string(gp::DIMENSIONS + 1));
        }

        std::complex<double> * data = static_cast<std::complex<double> *>(buf1.ptr);

        std::array<int,gp::DIMENSIONS> shape;


        for(int d=0;d<gp::DIMENSIONS;d++)
        {
            shape[d]=buf1.shape[d];

            //std::cout << shape[d] << std::endl;
            //std::cout << buf1.strides[d] << std::endl;
            
        }

        int pyComps=buf1.shape[gp::DIMENSIONS];

        if (pyComps != phi.nComp()/2)
        {
             throw std::runtime_error("Number of components expected " + 
             std::to_string(phi.nComp()/2));
        }

#if AMREX_SPACEDIM == 3    
            for (int k=region.minIndex(2);k<=region.maxIndex(2);k++)
	 		    for (int j=region.minIndex(1);j<=region.maxIndex(1);j++) 
	     	        for (int i=region.minIndex(0);i<=region.maxIndex(0);i++)
                        for(int c=0;c<pyComps;c+=1)
                        {
                            int index = i + j*shape[0] + k *shape[0] * shape[1] + c * shape[0] * shape[1] * shape[2];


                            phi_arr(i,j,k,2*c+1)=data[index].imag();
                            phi_arr(i,j,k,2*c)=data[index].real();
                        }
            
#endif




    } 

    auto filename = settings["folder"].get<std::string>() + std::string("/") + settings["name"].get<std::string>() ;

   amrex::VisMF::Write(phi, filename);

    //amrex::Finalize();
   }



std::vector<std::vector<std::complex<double> > > readMultifab( const json_t & settings   )
{
    initializer::getInstance().init();


    auto [baDomain , geom , dm, low_bc, high_bc] = gp::createGeometry(settings["geometry"]);

    auto phi = createMultifab(settings);

    gp::wavefunction wave(&phi,&geom);

    auto filename = settings["folder"].get<std::string>() + std::string("/") + settings["name"].get<std::string>() ;


    amrex::VisMF::Read(phi, filename);

    std::vector< std::vector<std::complex<double> > >data;
    
    for ( auto mfi = wave.beginAmrexIterator() ; mfi.isValid(); ++mfi )
    {
        // loop over box region
        auto region = wave[mfi];
        auto & phi_arr = region.getPhi<gp::array_t>();

        //std::cout <<"region " <<  region.size() << std::endl;

        std::vector<std::complex<double> > regionData(region.size()/2,std::complex<double>(0,0));


        std::array<int,gp::DIMENSIONS> shape {AMREX_D_DECL( region.size(0),region.size(1),region.size(2))};

        

#if AMREX_SPACEDIM == 3    
            for (int k=region.minIndex(2);k<=region.maxIndex(2);k++)
	 		    for (int j=region.minIndex(1);j<=region.maxIndex(1);j++) 
	     	        for (int i=region.minIndex(0);i<=region.maxIndex(0);i++)
                        for(int c=0;c<region.nComp()/2;c+=1)
                        {
                            int index = i + j*shape[0] + k *shape[0] * shape[1] + c * shape[0] * shape[1] * shape[2];

                            regionData[index]=std::complex<double>(phi_arr(i,j,k,2*c) ,phi_arr(i,j,k,2*c + 1) );
                        }            
#endif
        data.push_back(regionData);

    }
    

    

    
    return data;
    
};












