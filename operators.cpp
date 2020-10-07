#include "operators.h"
#include "stencils.h"
#include "gpExceptions.h"
#include "tools.h"
#include <AMReX_BCUtil.H>

op::op()
{
    for (int i=0;i<AMREX_SPACEDIM;i++)
    {        _bc_lo[i] = BC::PERIODIC;
             _bc_hi[i] = BC::PERIODIC;
    }
    _nComponents=1;

}

void op::define (Geometry & geom_ , BoxArray & ba_ , DistributionMapping & dm_, bc_t & bc_low, bc_t & bc_high)
{
    _geom=geom_;
    _ba = ba_;
    _dm= dm_;
    _bc_lo = bc_low;
    _bc_hi = bc_high;
    _bc_multifab= {toMultiFabBC(bc_low,bc_high)} ;

}

void amrexLaplacianOperator::define (Geometry & geom_ , BoxArray & ba_ , DistributionMapping & dm_, bc_t & bc_low, bc_t & bc_high)
{
    laplacianOperator::define( geom_, ba_ , dm_ , bc_low, bc_high) ;

    int coord = getGeometry().Coord();

    info.setMetricTerm(true);
    info.setMaxCoarseningLevel(0);

    linPoisson.define({getGeometry()}, {getBoxArray()}, {getDistributionMapping() }, info);

    //linPoissonReal.setNComp(nComp);

    std::array<LinOpBCType,AMREX_SPACEDIM> bc_lo;
    std::array<LinOpBCType,AMREX_SPACEDIM> bc_hi;

    for (int d=0;d<amrex::SpaceDim;d++)
        {
        std::cout << getLowBC()[0] << std::endl;
        bc_lo[d]=toLinOpBCType(getLowBC()[d]);
        bc_hi[d]=toLinOpBCType(getHiBC()[d]);
        }
    
    linPoisson.setDomainBC(bc_lo, bc_hi);
  
    
    linPoisson.setLevelBC(0, nullptr);

    linPoisson.setMaxOrder(order);

    if (lap != NULL)
    {
        delete lap;
    }
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

stencilLaplacianOperator::stencilLaplacianOperator(int order_) : 
        order(order_)
     {
         if (order != 2) {
             throw missingImplementation("Higher order laplacian not implemented");
         }
    }


template<int order>
 void stencilLaplacian<order>::define (Geometry & geom_ , BoxArray & ba_ , DistributionMapping & dmInit, bc_t & bc_low, bc_t &  bc_high)
{
    laplacianOperator::define(geom_,ba_,dmInit,bc_low,bc_high);
  
    BoxArray edgeBaX=  getBoxArray().surroundingNodes(0)  ;
    fluxx.define( edgeBaX, getDistributionMapping(), getNComp(), 0 );

    # if AMREX_SPACEDIM >= 2
    BoxArray edgeBaY=  getBoxArray().surroundingNodes(1)  ;
    fluxy.define( edgeBaY, getDistributionMapping(), getNComp(), 0 );
    #endif  

    # if AMREX_SPACEDIM >= 3
    BoxArray edgeBaZ=  getBoxArray().surroundingNodes(2)  ;
    fluxz.define( edgeBaZ, getDistributionMapping(), getNComp(), 0 );
    #endif 

}

template<int order>
 void stencilLaplacian<order>::apply(MultiFab & newMultiFab, MultiFab & oldMultiFab) 
 {
     auto & geom = getGeometry();
     if ( geom.Coord() == 2 )
     {
         applySpherical(newMultiFab,oldMultiFab);
     }
     else if (geom.Coord() == 0 )
     {
         applyCartesian(newMultiFab,oldMultiFab);
     }
     else
     {
         throw missingImplementation("Only cartesian and spherical laplacain stencils are supported");
     }
     
 }

template<int order>
 void stencilLaplacian<order>::applyCartesian(MultiFab & newMultiFab, MultiFab & oldMultiFab) 
 {
     auto & geom = getGeometry();
   // oldMultiFab.FillBoundary(getGeometry().periodicity());

    amrex::FillDomainBoundary(oldMultiFab, getGeometry(), getBCRec() );

	const Real* dx = geom.CellSize(); 
	const Real* prob_lo = geom.ProbLo();

    std::array<Real, AMREX_SPACEDIM> dxInverseSquare;

    for (int d=0;d<AMREX_SPACEDIM;d++)
    {
        dxInverseSquare[d]=1./(dx[d]*dx[d]);
    } 

   

	for ( MFIter mfi(newMultiFab); mfi.isValid(); ++mfi ) 
	{ 
	    const Box& bx = mfi.validbox(); 
	    const int* lo = bx.loVect(); 
	    const int *hi= bx.hiVect(); 
	    Array4< Real> const & newData = newMultiFab[mfi].array();
        Array4< Real> const & oldData = oldMultiFab[mfi].array();

#if    AMREX_SPACEDIM == 3
     for (int k=lo[2];k<=hi[2]  ;k++) 
        for (int j=lo[1];j<=hi[1]  ;j++) 
            for (int i=lo[0];i<=hi[0] ;i++) 
            {

                newData(i,j,k)= stencil1D<order,2>::x(oldData,i,j,k) *  dxInverseSquare[0] + stencil1D<order,2>::y(oldData,i,j,k) * dxInverseSquare[1] + stencil1D<order,2>::z(oldData,i,j,k) *  dxInverseSquare[2]   ;

            }
#endif

#if         AMREX_SPACEDIM == 1

            int j=0;
            int k=0;

            for (int i=lo[0];i<=hi[0] ;i++) 
            {
                newData(i,j,k)= stencil1D<order,2>::x(oldData,i,j,k) * dxInverseSquare[0];
                

            }
#endif

    }

 }


 template<int order>
 void stencilLaplacian<order>::applySpherical(MultiFab & newMultiFab, MultiFab & oldMultiFab) 
 {
     auto & geom = getGeometry();
     
    amrex::FillDomainBoundary(oldMultiFab, getGeometry(), getBCRec() );
    
	const Real* dx = geom.CellSize(); 
	const Real* prob_lo = geom.ProbLo();

    std::array<Real, AMREX_SPACEDIM> dxInverseSquare;

    for (int d=0;d<AMREX_SPACEDIM;d++)
    {
        dxInverseSquare[d]=1./(dx[d]*dx[d]);
    } 

	for ( MFIter mfi(newMultiFab); mfi.isValid(); ++mfi ) 
	{ 
	    const Box& bx = mfi.validbox(); 
	    const int* lo = bx.loVect(); 
	    const int *hi= bx.hiVect(); 
	    Array4< Real> const & newData = newMultiFab[mfi].array();
        Array4< Real> const & oldData = oldMultiFab[mfi].array();

#if         AMREX_SPACEDIM == 1

            int j=0;
            int k=0;

            for (int i=lo[0];i<=hi[0] ;i++) 
            {
                auto x = prob_lo[0] + (i+0.5)*dx[0];

                newData(i,j,k)= (stencil1D<order,2>::x(oldData,i,j,k) + 2*stencil1D<order,1>::x(oldData,i,j,k)/x *dx[0]  ) * dxInverseSquare[0];

            }
#endif

#if AMREX_SPACEDIM > 1
    throw invalidInput("Spherica dimensions only supported in dimensions 1.");
#endif



    }

 }


/*

template<int order>
 void stencilLaplacian<order>::apply2OrderFiniteVolume(MultiFab & newMultiFab, MultiFab & oldMultiFab) 
{   
    auto & geom = getGeometry();

	const Real* dx = geom.CellSize(); 
	const Real* prob_lo = geom.ProbLo(); 

    oldMultiFab.FillBoundary(geom.periodicity() );

	for ( MFIter mfi(newMultiFab); mfi.isValid(); ++mfi ) 
	{ 
	    const Box& bx = mfi.validbox(); 
	    const int* lo = bx.loVect(); 
	    const int *hi= bx.hiVect(); 
	    Array4< Real> const & newData = newMultiFab[mfi].array();
        Array4< Real> const & oldData = oldMultiFab[mfi].array();


        Array4<Real> const & fluxxData = fluxx[mfi].array();
        #if AMREX_SPACEDIM >=2
        Array4<Real> const & fluxyData = fluxy[mfi].array();
        #endif

        #if AREX_SPACEDIM>=3
        Array4<Real> const & fluxzData = fluxz[mfi].array();
        #endif

        // -------   Compute fluxes ----------------- s

#if AMREX_SPACEDIM == 1
		const int j=0;
		const int k=0;
        for (int i=lo[0];i<=hi[0] + 1;i++) 
			{
                fluxxData(i,j,k) = stencil1D<order,1>::x(oldData,i,j,k) / dx[0];

            }
                //-----------  evaluates the second derivatives -----s

        for (int i=lo[0];i<=hi[0] ;i++) 
			{
                newData(i,j,k)=   stencil1D<order,1>::x_nodal(fluxxData,i,j,k)/ dx[0]  
                ;
            }

#endif

#if AMREX_SPACEDIM == 3
     for (int k=lo[2];k<=hi[2]  ;k++) 
        for (int j=lo[1];j<=hi[1]  ;j++) 
            for (int i=lo[0];i<=hi[0] + 1;i++) 
  
    // flux in the x direction
    		{
                fluxxData(i,j,k) = stencil1D<order,1>::x(oldData,i,j,k) / dx[0];
            }
            
    // flux in the y direction
    for (int k=lo[2];k<=hi[2] ;k++)
     for (int j=lo[1];j<=hi[1] + 1 ;j++) 
        for (int i=lo[0];i<=hi[0] ;i++) 
			{
                fluxyData(i,j,k) = stencil1D<order,1>::y(oldData,i,j,k) / dx[1];
            }

    // flux in the z direction
    for (int k=lo[2];k<=hi[2] + 1 ;k++) 
     for (int j=lo[1];j<=hi[1]  ;j++) 
        for (int i=lo[0];i<=hi[0] ;i++) 
  			{
                fluxzData(i,j,k) = stencil1D<order,1>::z(oldData,i,j,k) / dx[2];
            }

        for (int k=lo[2];k<=hi[2];k++) 
            for (int j=lo[1];j<=hi[1];j++)
                for (int i=lo[0];i<=hi[0] ;i++) 
			    {
                    newData(i,j,k)=   stencil1D<order,1>::x_nodal(fluxxData,i,j,k)/ dx[0]  +  stencil1D<order,1>::y_nodal(fluxyData,i,j,k)/ dx[1] + 
                     stencil1D<order,1>::z_nodal(fluxzData,i,j,k)/ dx[2]
                ;
                }
#endif

    }


}


*/



template class stencilLaplacian<2>;
