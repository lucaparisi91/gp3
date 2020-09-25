#ifndef TOOLS_H
#define TOOLS_H
#include <AMReX_PlotFileUtil.H>
#include <nlohmann/json.hpp>
using json_t = nlohmann::json ;
using namespace amrex;



#define LOOP3D( state , geom  ) \
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

#define ENDLOOP3D \
			 }\
	} \
	} 





Real norm( const MultiFab & phi_real , const MultiFab & phi_imag,  const Geometry & geom, int component=0);

std::tuple< BoxArray , Geometry , DistributionMapping   >  
createGeometry( const json_t & settings);

template<class evaluator_t >
void fill( MultiFab & state, Geometry & geom,  const evaluator_t & evaluator )
{

  LOOP3D( state, geom  )

   auto x = prob_lo[0] +  (i + 0.5) * dx[0] ;
	 auto y= prob_lo[1] + (j + 0.5) * dx[1];
	 auto z= prob_lo[2] + (k + 0.5) * dx[2];		

   data(i,j,k) = evaluator(x,y,z);

  ENDLOOP3D

}

class initializer
{
    bool isInitialized;
    static initializer *s_instance;
    initializer() : isInitialized(false ) {}
  

  public:
    void init()
    {
        if (!isInitialized)
        {
            amrex::Initialize(MPI_COMM_WORLD);
            isInitialized = true;
        }
    }



    static initializer  & instance()
    {
        if (!s_instance)
          s_instance = new initializer;
        return *s_instance;
    }
};

// Allocating and initializing GlobalClass's
// static data member.  The pointer is being
// allocated - not the object inself.

#endif