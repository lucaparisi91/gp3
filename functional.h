#include "evaluate.h"
#include <AMReX_PlotFileUtil.H>
#include "operators.h"






#define EVALUATION_LOOP3D( state_new_real, state_new_imag, state_old_real, state_old_imag ,  geom  ) { \
	const Real* dx = geom.CellSize(); \
	const Real* prob_lo = geom.ProbLo(); \
	for ( MFIter mfi(state_new_real); mfi.isValid(); ++mfi ) \
	{ \
	    const Box& bx = mfi.validbox(); \
	    const int* lo = bx.loVect(); \
	    const int *hi= bx.hiVect(); \
		Array4< const Real> const & phi_old_real = state_old_real[mfi].const_array(); \
	    Array4< Real> const & phi_new_real = state_new_real[mfi].array(); \
	    Array4< const Real> const & phi_old_imag =  state_old_imag[mfi].const_array();\
	    Array4< Real> const & phi_new_imag = state_new_imag[mfi].array(); \
		for (int k=lo[2];k<=hi[2];k++) \
	 		for (int j=lo[1];j<=hi[1];j++) \
	     		for (int i=lo[0];i<=hi[0];i++) \
	     	{ 



#define END_EVALUATION_LOOP3D }\
			 }\
			 }











class functional
{
public:
    functional() : lap(NULL) {}
    virtual void define( Geometry & geom_ , BoxArray & box_, DistributionMapping & dm_ );

    virtual void evaluate(
	MultiFab & state_new_real, MultiFab & state_new_imag,
	MultiFab & state_old_real,  MultiFab & state_old_imag , Real time=0) = 0;



    auto & getGeometry() {return _geom;} 
    auto & getBoxArray() {return _box;}
    auto & getDistributionMapping() {return _dm;}

    auto & getLaplacian() {return *lap;}

    ~functional();

    private:
    Geometry _geom;
    BoxArray _box;
    DistributionMapping _dm;
    int _laplacianOrder;
    laplacianOperator *lap;
};

class harmonicFunctional : public functional
{
    public:

    harmonicFunctional(Real omega = 1) : prefactor(0.5 * omega*omega) {};

	harmonicFunctional(const json_t & j ) : harmonicFunctional( j["omega"].get<Real>() ) {}


    virtual void evaluate(
	MultiFab & state_new_real, MultiFab & state_new_imag,
	MultiFab & state_old_real,  MultiFab & state_old_imag, Real time =0 ) ;

	static std::string name() {return "harmonic";}
	
    private:

    Real prefactor;

};




