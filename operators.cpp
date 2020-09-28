#include "operators.h"

void laplacianOperator::define (Geometry & geom_ , BoxArray & ba_ , DistributionMapping & dm_)
{
    _geom=geom_;
    _ba = ba_;
    _dm= dm_;
}

void amrexLaplacianOperator::define (Geometry & geom_ , BoxArray & ba_ , DistributionMapping & dm_)

{
    laplacianOperator::define( geom_, ba_ , dm_ ) ;


    Vector<BCRec> bc(Ncomp);

    for (int n = 0; n < Ncomp; ++n)
        {
            for (int d=0;d<amrex::SpaceDim ; d++)
            {    
                bc[n].setLo(d, BCType::int_dir);  
                bc[n].setHi(d, BCType::int_dir);
        }
    }

    //info.setMetricTerm(true);
    info.setMaxCoarseningLevel(0);

    linPoisson.define({getGeometry()}, {getBoxArray()}, {getDistributionMapping() }, info);

    //linPoissonReal.setNComp(nComp);
    
    linPoisson.setMaxOrder(order);

    std::array<LinOpBCType,AMREX_SPACEDIM> bc_lo;
    std::array<LinOpBCType,AMREX_SPACEDIM> bc_hi;


    for (int d=0;d<amrex::SpaceDim;d++)
        {
        bc_lo[d]=LinOpBCType::Periodic;
        bc_hi[d]=LinOpBCType::Periodic;
        }
    
    linPoisson.setDomainBC(bc_lo, bc_hi);
  
    linPoisson.setLevelBC(0,nullptr);

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
