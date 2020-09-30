#include "operators.h"
#include "stencils.h"
#include "gpExceptions.h"


void laplacianOperator::define (Geometry & geom_ , BoxArray & ba_ , DistributionMapping & dm_)
{
    _geom=geom_;
    _ba = ba_;
    _dm= dm_;
}

void amrexLaplacianOperator::define (Geometry & geom_ , BoxArray & ba_ , DistributionMapping & dm_)

{
    laplacianOperator::define( geom_, ba_ , dm_ ) ;

    int coord = getGeometry().Coord();


    Vector<BCRec> bc(Ncomp);

    for (int n = 0; n < Ncomp; ++n)
        {
            for (int d=0;d<amrex::SpaceDim ; d++)
            {    
                bc[n].setLo(d, BCType::int_dir);  
                bc[n].setHi(d, BCType::int_dir);
            }

            if (coord == 2)
            {   
                bc[n].setLo(0,BCType::foextrap);
                bc[n].setHi(0,BCType::foextrap);
            }

    }

    info.setMetricTerm(true);
    info.setMaxCoarseningLevel(0);

    linPoisson.define({getGeometry()}, {getBoxArray()}, {getDistributionMapping() }, info);

    //linPoissonReal.setNComp(nComp);
    

    std::array<LinOpBCType,AMREX_SPACEDIM> bc_lo;
    std::array<LinOpBCType,AMREX_SPACEDIM> bc_hi;

    for (int d=0;d<amrex::SpaceDim;d++)
        {
        bc_lo[d]=LinOpBCType::Periodic;
        bc_hi[d]=LinOpBCType::Periodic;
        }
    
    if (coord == 2)
    {
        bc_lo[0]=LinOpBCType::Neumann;
        bc_hi[0]=LinOpBCType::Neumann;
    }
    
    linPoisson.setDomainBC(bc_lo, bc_hi);
  
    linPoisson.setLevelBC(0,nullptr);

    linPoisson.setMaxOrder(4);
    
    lap= new MLMG(linPoisson);

}

amrexLaplacianOperator::~amrexLaplacianOperator()
{
    if (lap != NULL)
    {
        delete lap;
    }
    

}


void amrexLaplacianOperator::apply(MultiFab & newMultiFab, MultiFab & oldMultiFab)
{
    
    lap->apply({ &newMultiFab},{ &oldMultiFab} ) ;
}

void stencilLaplacianOperator::apply2OrderSpherical(MultiFab & state, MultiFab & stateOld )
{
#if AMREX_SPACEDIM == 1

    auto & geom = getGeometry();

    const Real* dx = geom.CellSize();
	const Real* prob_lo = geom.ProbLo();
    const int *  hiDomain = geom.Domain().hiVect();

	for ( MFIter mfi(state); mfi.isValid(); ++mfi ) 
	{ 
	    const Box& bx = mfi.validbox();
	    const int* lo = bx.loVect(); 
	    const int *hi= bx.hiVect(); 
	    Array4< Real> const & stateFab = state[mfi].array();
		Array4< Real> const & stateOldFab = stateOld[mfi].array();
		const int j=0;
		const int k=0;
        Real x=0;

        // loop over the inner part of the box
        int iLow=std::max( lo[0] , order - 1 );
        int iHi=std::min( hi[0] , hiDomain[0] - order + 1 );
        

		for (int i=iLow;i<=iHi;i++) 
			{
                x = prob_lo[0] +  (i + 0.5) * dx[0];

                stateFab(i,j,k)=laplacianSphericalSymm<2,CENTRAL>::call(stateOldFab,i,j,k,dx,x);
            }
        
        // loop over the left boundary
        for(int i=lo[0];i<iLow;i++)
        {
            x = prob_lo[0] +  (i + 0.5) * dx[0];
            stateFab(i,j,k)=laplacianSphericalSymm<2,FORWARD>::call(stateOldFab,i,j,k,dx,x);
        }

        // loop over the right boundary
        for(int i=iHi + 1 ;i<=std::min(hiDomain[0],hi[0]);i++)
        {
            x = prob_lo[0] +  (i + 0.5) * dx[0];

            stateFab(i,j,k)=laplacianSphericalSymm<2,FORWARD>::call(stateOldFab,i,j,k,dx,x);
        }
    }

#endif





}

void stencilLaplacianOperator::apply(MultiFab & state, MultiFab & stateOld )
{
    if (getGeometry().Coord() == 2 and (AMREX_SPACEDIM == 1 ) )
    {
        apply2OrderSpherical(state,stateOld);
    }
    else
    {
        throw missingImplementation("Non spherical higher order stencil laplacian not implemented yet");

    }

}