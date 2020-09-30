#ifndef OPERATORS_H
#define OPERATORS_H


#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>
#include <AMReX_PlotFileUtil.H>
#include "traits.h"

using namespace amrex;


class op
{
    public:
    op( ) {}
    virtual void define (Geometry & geom_ , BoxArray & ba_ , DistributionMapping & dmInit);

    virtual void apply(MultiFab & newMultiFab, MultiFab & oldMultiFab) = 0;

    auto & getGeometry() {return _geom;}
    auto & getBoxArray() {return _ba;}
    auto & getDistributionMapping() {return _dm;}

    private:
    Geometry _geom;
    BoxArray _ba;
    DistributionMapping _dm;
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
        virtual void define (Geometry & geom_ , BoxArray & ba_ , DistributionMapping & dmInit);
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

#endif