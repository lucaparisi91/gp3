#ifndef OPERATORS_H
#define OPERATORS_H


#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>
#include <AMReX_PlotFileUtil.H>

using namespace amrex;


class laplacianOperator
{
    public:
    laplacianOperator( ) {}
    virtual void define (Geometry & geom_ , BoxArray & ba_ , DistributionMapping & dmInit);

    virtual void apply(MultiFab & newMultiFab, MultiFab & oldMultiFab) = 0;

    auto & getGeometry() {return _geom;}
    auto & getBoxArray() {return _ba;}
    auto & getDistributionMapping() {return _dm;}

    ~laplacianOperator() {}


    private:
    Geometry _geom;
    BoxArray _ba;
    DistributionMapping _dm;
};

class amrexLaplacianOperator : public laplacianOperator
{
    public:     
        amrexLaplacianOperator() : order(2) ,Ncomp(1), lap(NULL ) , laplacianOperator() {}
        virtual void define (Geometry & geom_ , BoxArray & ba_ , DistributionMapping & dmInit);
        virtual void apply(MultiFab & newMultiFab, MultiFab & oldMultiFab) ;

         ~amrexLaplacianOperator();
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
    stencilLaplacianOperator(int order_ = 2) : order(order_), laplacianOperator::laplacianOperator(){} ;

    virtual void apply(MultiFab & newMultiFab, MultiFab & oldMultiFab) ;
    void apply2OrderSpherical(MultiFab & newMultiFab, MultiFab & oldMultiFab) ;


    private:

    int order;
};

#endif