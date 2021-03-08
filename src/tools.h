#include "traits.h"

namespace gp
{


amrex::BCRec toMultiFabBC(const std::array<BC,AMREX_SPACEDIM> & bc_low, const std::array<BC,AMREX_SPACEDIM> & bc_high  );

std::tuple<std::array<BC,AMREX_SPACEDIM>, std::array<BC,AMREX_SPACEDIM> > readBC(const json_t & j);


};
