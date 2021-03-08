#ifndef OPERATORS_H
#define OPERATORS_H

#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>
#include <AMReX_PlotFileUtil.H>
#include "src/traits.h"
#include "tools.h"

using namespace amrex;

class op
{
    public:
    op( ) ;
    virtual void define (Geometry & geom_ , BoxArray & ba_ , DistributionMapping & dmInit, bc_t & bc_low, bc_t &  bc_high);

    virtual void apply(MultiFab & newMultiFab, MultiFab & oldMultiFab) = 0;

    auto & getGeometry() {return _geom;}
    auto & getBoxArray() {return _ba;}
    auto & getDistributionMapping() {return _dm;}

    const auto & getLowBC() const  {return _bc_lo;}
    const auto & getHiBC() const {return _bc_hi;}

    auto getNComp() const {return _nComponents;}

    auto &  getBCRec() {return _bc_multifab;}
    private:
    Geometry _geom;
    BoxArray _ba;
    DistributionMapping _dm;
    int _nComponents;

    std::array<BC, AMREX_SPACEDIM> _bc_lo;
    std::array<BC, AMREX_SPACEDIM> _bc_hi;

    Vector< BCRec > _bc_multifab;

};

class laplacianOperator : public op
{
public:

    laplacianOperator() : op::op() {}

};

class amrexLaplacianOperator : public laplacianOperator
{
    public:    
        amrexLaplacianOperator(const json_t & j ) : 
        amrexLaplacianOperator::amrexLaplacianOperator(j["order"].get<int>() ) {
        }

        amrexLaplacianOperator(int order_=2) : order(order_) ,Ncomp(1), lap(NULL ) , laplacianOperator() {}
        virtual void define (Geometry & geom_ , BoxArray & ba_ , DistributionMapping & dmInit, bc_t & bc_low, bc_t & bc_high);


        virtual void apply(MultiFab & newMultiFab, MultiFab & oldMultiFab) ;
         ~amrexLaplacianOperator();

         
    

         static std::string name() {return "amrexLaplacian";}
    private :
        LPInfo info;
        int order ;
        MLPoisson linPoisson;
        int Ncomp ;
        MLMG *lap;
};

class stencilLaplacianOperator : public laplacianOperator
{
    public:

    stencilLaplacianOperator(const json_t & j ) : 
        stencilLaplacianOperator::stencilLaplacianOperator(j["order"].get<int>() ) {
        }

    stencilLaplacianOperator(int order_ = 2) ;


    virtual void apply(MultiFab & newMultiFab, MultiFab & oldMultiFab) ;
    void apply2OrderSpherical(MultiFab & newMultiFab, MultiFab & oldMultiFab) ;

    static std::string name() {return "stencilLaplacian";}

    private:

    int order;
};




template<int order>
class stencilLaplacian : public laplacianOperator
{
    public:
    stencilLaplacian() : laplacianOperator::laplacianOperator() {}

    stencilLaplacian(const json_t & j) : stencilLaplacian() {}

    virtual void define (Geometry & geom_ , BoxArray & ba_ , DistributionMapping & dmInit, bc_t & bc_low, bc_t &  bc_high);

    virtual void apply(MultiFab & newMultiFab, MultiFab & oldMultiFab) ;

    virtual void applySpherical(MultiFab & newMultiFab, MultiFab & oldMultiFab) ;
    virtual void applyCartesian(MultiFab & newMultiFab, MultiFab & oldMultiFab) ;

    auto static name(){ return  "stencilLaplacian" + std::to_string(order ) ;}

    private:
    /* temporary storage for fluxes*/
    MultiFab fluxx ;
    MultiFab fluxy ;
    MultiFab fluxz ;
};


#endif