#include "tools.h"
#include "gpExceptions.h"

namespace gp
{
    using namespace amrex;


std::tuple<std::array<BC,AMREX_SPACEDIM>, std::array<BC,AMREX_SPACEDIM> > readBC( const json_t & j)
{
    std::array<BC,AMREX_SPACEDIM> low_bc;
    std::array<BC,AMREX_SPACEDIM> high_bc;

    for (int d=0;d<AMREX_SPACEDIM;d++)
    {
        BC bc = BC::PERIODIC;
        auto bc_type = j[d].get<std::string>(); 
        if ( bc_type == "neumann" )
        {
            bc=BC::NEUMANNN;
        }
        else
        { if (bc_type == "drichlet")
            {
                bc=BC::DRICHLET;
            }
            else if (bc_type != "periodic")
            {
                throw invalidInput("Unkown bc " + bc_type);
            }
        }
        low_bc[d]=bc;
        high_bc[d] = bc;
    }

    return std::make_tuple(low_bc,high_bc);
}

auto toMultiFabBC(const BC & bc)
{
    if (bc== BC::PERIODIC)
    {
        return BCType::int_dir;
    }
    else if (bc == BC::NEUMANNN)
    {
        return BCType::reflect_even;
    }
    else if (bc == BC::DRICHLET)
    {
        return BCType::ext_dir;
    }
    else
    {
        throw invalidInput("Unkown boundary condition");
    }

}


BCRec toMultiFabBC(const std::array<BC,AMREX_SPACEDIM> & bc_low, const std::array<BC,AMREX_SPACEDIM> & bc_high  )
{
  BCRec bc;
    
    for (int d=0;d<AMREX_SPACEDIM;d++)
    {  
        bc.setLo(d, toMultiFabBC(bc_low[d]) ); 
        bc.setHi(d, toMultiFabBC(bc_high[d]));
    }

    return bc;
}



};