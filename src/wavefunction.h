#include "geometry.h"

namespace gp
{

    using box=amrex::Box;
    using geometry=amrex::Geometry;
    
    using multiFab = amrex::MultiFab;
    using fab = amrex::FArrayBox;
    using array_t = amrex::Array4<Real> ;
    using const_array_t = amrex::Array4<const Real>; 
    using smallVector = std::array<int,DIMENSIONS>;



    struct wavefunctionRegion;


    struct wavefunction
    {
        // contains the wavefunction at a single level of refinement
        wavefunction(multiFab * phiNew_, multiFab * phiOld_, const geometry * geom_) : phiNew(phiNew_),phiOld(phiOld_),geom(geom_)   {
        }

        const auto & getGeometry() const {return *geom;};

        const auto & getPhiNew() const {return *phiNew;}
        auto & getPhiNew() {return *phiNew;}

        auto & getPhiOld() {return *phiOld;}
        const auto & getPhiOld() const {return *phiOld;}

        wavefunctionRegion operator[](const amrex::MFIter & it);

        amrex::MFIter beginAmrexIterator() {return amrex::MFIter(getPhiNew() ) ;} 

        void fillBoundaries() // fills all boundary. Only pbc supported at the moment
        {
            phiNew->FillBoundary(geom->periodicity());
        }


        void swapOldAndNew(){
            std::swap(phiNew,phiOld);
        }


        private:

        multiFab * phiNew;
        multiFab * phiOld;
        const geometry * geom;

    };

    struct  wavefunctionRegion
    {
        wavefunctionRegion( box boxRegion_,fab * phiNew_, fab* phiOld_,const geometry* geom_ ) :
        boxRegion(boxRegion_),

        phiNew(phiNew_),
        phiOld(phiOld_),
        geom(geom_)
        
         {
             phiNewArr=phiNew->array();
             phiOldArr=phiOld->array();
             
             
         }

        int minIndex(int d) {return boxRegion.loVect()[d]; }
        int maxIndex(int d){return boxRegion.hiVect()[d];}

        smallVector minIndex();
        smallVector maxIndex();


        const auto & getBox(){return boxRegion;}
        
        auto & getPhiNewAsFab(){return *phiNew;}
        const auto & getPhiNewAsFab() const {return *phiNew;}
        auto & getPhiOldAsFab(){return *phiOld;}
        const auto & getPhiOldAsFab() const {return *phiOld;}


        // enable get access to wavefunction data
        template<class T, std::enable_if_t<std::is_same<T,fab>::value > * = nullptr   >
        auto & getPhiNew() { return getPhiNewAsFab();}

        template<class T, std::enable_if_t<std::is_same<T,array_t>::value > * = nullptr  >
        const auto & getPhiNew() { return phiNewArr;}


        template<class T, std::enable_if_t<std::is_same<T,array_t>::value > * = nullptr  >
        const auto & getPhiOld() const { return phiOldArr;}

        const auto & getGeometry() const {return *geom;};

    private:
        box boxRegion;
        fab * phiNew;
        fab * phiOld;
        const geometry * geom;

        gp::array_t phiNewArr;
        gp::array_t phiOldArr;        

    };

};

