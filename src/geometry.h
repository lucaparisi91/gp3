
#include "gpExceptions.h"
#include "traits.h"


namespace gp
{

std::tuple<  amrex::Geometry  , std::array<BC,AMREX_SPACEDIM>,   std::array<BC,AMREX_SPACEDIM>  >
createGeometry( const json_t & settings);
amrex::Box createDomainBox(const json_t & settings);
amrex::BoxArray createGrids(const json_t & j);
amrex::BoxArray createBoxes(const json_t & j);


json_t toJson(const amrex::Geometry & geo);
json_t toJson(const amrex::BoxArray & box,const amrex::DistributionMapping & dm);
json_t toJson(const amrex::BoxArray & box);




};

