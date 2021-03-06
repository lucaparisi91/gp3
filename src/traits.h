#ifndef TRAITS
#define TRAITS

#include <nlohmann/json.hpp>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>


namespace gp
{
enum BC { PERIODIC = 0, DRICHLET = 1 , NEUMANNN = 2 };

using Real=double;
using json_t = nlohmann::json ;
constexpr int DIMENSIONS = AMREX_SPACEDIM;

}



#endif