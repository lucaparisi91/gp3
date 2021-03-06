
#include "gpExceptions.h"
#include "traits.h"


namespace gp
{





std::tuple< amrex::BoxArray , amrex::Geometry , amrex::DistributionMapping  , std::array<BC,AMREX_SPACEDIM>,   std::array<BC,AMREX_SPACEDIM>  >

createGeometry( const json_t & settings);
};
