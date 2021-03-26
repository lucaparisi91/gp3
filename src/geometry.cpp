

#include "gpExceptions.h"
#include "geometry.h"
#include "tools.h"
#include <AMReX_VisMF.H>


namespace gp
{


amrex::BoxArray createBoxes(const json_t & j)
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
};


amrex::Box createDomainBox(const json_t & settings)
{
    std::array<size_t,AMREX_SPACEDIM> shape;

    for (int i=0;i<AMREX_SPACEDIM;i++)
    {
        shape[i] = settings["shape"][i].get<int>();

    }

    amrex::IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
    amrex::IntVect dom_hi(AMREX_D_DECL( shape[0]-1, shape[1]-1, shape[2]-1));


    amrex::Box domain(dom_lo, dom_hi);

    return domain;


}


amrex::BoxArray createGrids(const json_t & j)
{
    if (j.find("boxes") != j.end() )
    {
        return createBoxes(j["boxes"]);
    }
    else
    {
        amrex::BoxArray ba;
        auto domain=createDomainBox(j["geometry"]);
        ba.define(domain);
        if (j.contains("maxGridSize") )
        {
            int max_grid_size = j["maxGridSize"].get<int>();
            ba.maxSize(max_grid_size);
        }
        return ba;
    }

};





std::tuple< amrex::Geometry   , std::array<BC,AMREX_SPACEDIM>,   std::array<BC,AMREX_SPACEDIM>  >
createGeometry( const json_t & settings)
{
    amrex::BoxArray ba;
    amrex::Geometry geom;
    amrex::Vector<int> is_periodic(AMREX_SPACEDIM,1);

    std::array<Real,AMREX_SPACEDIM> lower_edges;
    std::array<Real,AMREX_SPACEDIM> higher_edges;
    
    std::string coordinates = settings["coordinates"].get<std::string>();


    for (int i=0;i<AMREX_SPACEDIM;i++)
    {
        lower_edges[i]=settings["domain"][i][0].get<Real>();
        higher_edges[i]=settings["domain"][i][1].get<Real>();

    }
    
    auto domain=createDomainBox(settings);



    amrex::RealBox real_box({AMREX_D_DECL( lower_edges[0],lower_edges[1],lower_edges[2]) },
                         {AMREX_D_DECL( higher_edges[0], higher_edges[1], higher_edges[2] )});
                        
    int coord = 0;

    if (coordinates == "cartesian")
    {
        coord=0;
    }
    else if (coordinates == "spherical")
    {
        coord=2;

        is_periodic[0]=0;

    }
    else
    {
        throw missingImplementation("Unkown coordinates :" + coordinates );
    }

    std::array<BC,AMREX_SPACEDIM> low_bc , high_bc ;

    if (settings.contains("bc") )
    {
        std::tie(low_bc, high_bc) = readBC(settings["bc"]);

        for (int d=0;d<AMREX_SPACEDIM;d++)
        {
            if ( settings["bc"][d] != "periodic" )
            {
                is_periodic[d]=0;
            }    
        }

    }
    else
    {
        auto bc_settings =
        #if AMREX_SPACEDIM == 1     
        R"( ["periodic"   ])"_json
        #endif
         #if AMREX_SPACEDIM == 3     
        R"( ["periodic" , "periodic" , "periodic"  ])"_json
        #endif
        ;
        std::tie(low_bc, high_bc) = readBC(bc_settings);

    }
    

    geom.define(domain,&real_box,coord,is_periodic.data());

    
    return { geom, low_bc, high_bc} ;

}


}