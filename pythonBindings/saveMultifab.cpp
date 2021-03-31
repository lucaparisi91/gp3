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




void saveMultifab( py::list initialConditions , const json_t & settings   )
{
    initializer::getInstance().init();

    auto [ geom , low_bc, high_bc] = gp::createGeometry(settings["geometry"]);

    auto phi = gp::createMultiFab(settings);

    std::string name=settings["name"].get<std::string>();

    gp::wavefunction wave(&phi,&geom,name);
    
    auto nGhosts = settings["nGhosts"].get<std::vector<int> >();




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
            for (int k=region.minIndex(2)-nGhosts[2];k<=region.maxIndex(2) + nGhosts[2];k++)
	 		    for (int j=region.minIndex(1)-nGhosts[1];j<=region.maxIndex(1)+nGhosts[1];j++) 
	     	        for (int i=region.minIndex(0)-nGhosts[0];i<=region.maxIndex(0)+nGhosts[0];i++)
                        for(int c=0;c<pyComps;c+=1)
                        {
                            int ii=i - region.minIndex(0) + nGhosts[0];
                            int jj=j - region.minIndex(1) + nGhosts[1];
                            int kk=k - region.minIndex(2) + nGhosts[2];

                            int index = ii + jj*shape[0] + kk *shape[0] * shape[1] + c * shape[0] * shape[1] * shape[2];




                            phi_arr(i,j,k,2*c+1)=data[index].imag();
                            phi_arr(i,j,k,2*c)=data[index].real();
                        }
            
#endif




    } 

    //auto filename = settings["folder"].get<std::string>() + std::string("/") + settings["name"].get<std::string>() ;

    // amrex::VisMF::Write(phi, filename);
    std::string folder = settings["folder"].get<std::string>();
    wave.save(folder);
    //amrex::Finalize();
   }



std::vector<std::vector<std::complex<double> > > readMultifab( const json_t & settings   )
{
    initializer::getInstance().init();


    auto [ geom , low_bc, high_bc] = gp::createGeometry(settings["geometry"]);


    auto phi = gp::createMultiFab(settings);

    gp::wavefunction wave(&phi,&geom);

    auto nGhosts = settings["nGhosts"].get<std::vector<int> >();

    

    std::vector< std::vector<std::complex<double> > >data;
    
    for ( auto mfi = wave.beginAmrexIterator() ; mfi.isValid(); ++mfi )
    {
        // loop over box region
        auto region = wave[mfi];
        auto & phi_arr = region.getPhi<gp::array_t>();

        //std::cout <<"region " <<  region.size() << std::endl;
        
        size_t grownSize=1;

        std::array<int,gp::DIMENSIONS> shape {AMREX_D_DECL( region.size(0) + 2*nGhosts[0],region.size(1)+2*nGhosts[1],region.size(2) + 2*nGhosts[2])};


        for (int d=0;d<gp::DIMENSIONS;d++)
        {
            grownSize*=shape[d];
        }        


        grownSize*=phi.nComp()/2;

        std::vector<std::complex<double> > regionData(grownSize,std::complex<double>(0,0));

        std::cout << "reading..." << std::endl;


#if AMREX_SPACEDIM == 3    
    for (int k=region.minIndex(2)-nGhosts[2];k<=region.maxIndex(2) + nGhosts[2];k++)
	 		    for (int j=region.minIndex(1)-nGhosts[1];j<=region.maxIndex(1)+nGhosts[1];j++) 
	     	        for (int i=region.minIndex(0)-nGhosts[0];i<=region.maxIndex(0)+nGhosts[0];i++)
                        for(int c=0;c<phi.nComp()/2;c+=1)
                        {
                            int ii=i - region.minIndex(0) + nGhosts[0];
                            int jj=j - region.minIndex(1) + nGhosts[1];
                            int kk=k - region.minIndex(2) + nGhosts[2];

                            int index = ii + jj*shape[0] + kk *shape[0] * shape[1] + c * shape[0] * shape[1] * shape[2];



                            //std::cout << phi_arr(i,j,k,2*c) << std::endl; 
                            regionData[index]=std::complex<double>(phi_arr(i,j,k,2*c) ,phi_arr(i,j,k,2*c + 1) );


                        }            
#endif
        data.push_back(regionData);

    }
    

    
    
    
    return data;
    
};












