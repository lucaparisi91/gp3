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
    
    linPoisson.setMaxOrder(2);

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
