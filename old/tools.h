#ifndef TOOLS_H
#define TOOLS_H
#include <AMReX_PlotFileUtil.H>
#include <nlohmann/json.hpp>
using json_t = nlohmann::json ;
using namespace amrex;
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
namespace py = pybind11;
#include <AMReX_MLMG.H>
#include <AMReX_BCRec.H>



enum BC { PERIODIC = 0, DRICHLET = 1 , NEUMANNN = 2 };

LinOpBCType toLinOpBCType( const BC & bc );

using bc_t = std::array<BC,AMREX_SPACEDIM> ;

std::tuple<std::array<BC,AMREX_SPACEDIM>, std::array<BC,AMREX_SPACEDIM> > readBC(const json_t & j);

#if AMREX_SPACEDIM == 3
#define LOOP( state , geom  ) \
	{ \
	const Real* dx = geom.CellSize(); \
	const Real* prob_lo = geom.ProbLo(); \
	for ( MFIter mfi(state); mfi.isValid(); ++mfi ) \
	{ \
	    const Box& bx = mfi.validbox(); \
	    const int* lo = bx.loVect(); \
	    const int *hi= bx.hiVect(); \
	    Array4< Real> const & data = state[mfi].array();\
		for (int k=lo[2];k<=hi[2];k++) \
	 		for (int j=lo[1];j<=hi[1];j++) \
	     		for (int i=lo[0];i<=hi[0];i++) \
	     	{ \

#define ENDLOOP \
			 }\
	} \
	} 


#endif


#if AMREX_SPACEDIM==1

#define LOOP( state , geom  ) \
	{ \
	const Real* dx = geom.CellSize(); \
	const Real* prob_lo = geom.ProbLo(); \
	for ( MFIter mfi(state); mfi.isValid(); ++mfi ) \
	{ \
	    const Box& bx = mfi.validbox(); \
	    const int* lo = bx.loVect(); \
	    const int *hi= bx.hiVect(); \
	    Array4< Real> const & data = state[mfi].array();\
		const int j=0;\
		const int k=0;\
		for (int i=lo[0];i<=hi[0];i++) \
			{

#define LOOP_2ARG( state , state2,  geom  ) \
	{ \
	const Real* dx = geom.CellSize(); \
	const Real* prob_lo = geom.ProbLo(); \
	for ( MFIter mfi(state); mfi.isValid(); ++mfi ) \
	{ \
	    const Box& bx = mfi.validbox(); \
	    const int* lo = bx.loVect(); \
	    const int *hi= bx.hiVect(); \
	    Array4< Real> const & data = state[mfi].array();\
		Array4< Real> const & data2 = state[mfi].array();\
		const int j=0;\
		const int k=0;\
		for (int i=lo[0];i<=hi[0];i++) \
			{

#define ENDLOOP \
			 }\
	} \
	} 

#endif

Real norm( const MultiFab & phi_real , const MultiFab & phi_imag,  const Geometry & geom, int component=0);



void normalize(  MultiFab & phi_real ,  MultiFab & phi_imag,  const Geometry & geom, Real N=1);



std::tuple< BoxArray , Geometry , DistributionMapping  , std::array<BC,AMREX_SPACEDIM> , std::array<BC,AMREX_SPACEDIM> >  
createGeometry( const json_t & settings);

template<class evaluator_t >
void fill( MultiFab & state, Geometry & geom,  const evaluator_t & evaluator )
{
  LOOP( state, geom  )

#if AMREX_SPACEDIM == 3
   auto x = prob_lo[0] +  (i + 0.5) * dx[0] ;
   auto y= prob_lo[1] + (j + 0.5) * dx[1];
   auto z= prob_lo[2] + (k + 0.5) * dx[2];		

   data(i,j,k) = evaluator(x,y,z);
#endif

#if AMREX_SPACEDIM == 1
   auto x = prob_lo[0] +  (i + 0.5) * dx[0] ;
   data(i,j,k,0) = evaluator(x);
#endif



  ENDLOOP


}

void fill(MultiFab & realState, MultiFab & imagState, py::array_t<std::complex<Real> > initialCondition , Geometry & geom);

BCRec toMultiFabBC(const std::array<BC,AMREX_SPACEDIM> & bc_low, const std::array<BC,AMREX_SPACEDIM> & bc_high  );

#endif