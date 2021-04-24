#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include "geometry.h"


namespace gp
{
    using box=amrex::Box;
    using geometry=amrex::Geometry;
    using multiFab = amrex::MultiFab;
    using fab = amrex::FArrayBox;
    using array_t = amrex::Array4<Real> ;
    using constArray_t = amrex::Array4< const Real> ;

    using const_array_t = amrex::Array4<const Real>; 
    using smallVector = std::array<int,DIMENSIONS>;


    struct wavefunctionRegion;
    struct constWavefunctionRegion;

    amrex::MultiFab  createMultiFab(const json_t & settings);

    
    struct wavefunction
    {
        // contains the wavefunction at a single level of refinement
        wavefunction(multiFab * phi, const geometry * geom_, std::string name="phi") : _phi(phi),geom(geom_),_isReal(false),_name(name)   {

            if ( not isReal() )
            {
                assert( nComp()%2 == 0 );
            } 

        }

        const auto & getGeometry() const {return *geom;};

        auto & getPhi() {return *_phi;}
        
        const auto & getPhi() const {return *_phi;}

        const auto nGhosts() const {return _phi->nGrow();}


        wavefunctionRegion operator[](const amrex::MFIter & it) ;



        constWavefunctionRegion operator[](const amrex::MFIter & it) const ;

        amrex::MFIter beginAmrexIterator() {return amrex::MFIter(getPhi() ) ;}

        void fillBoundaries() // fills all boundary. Only pbc supported at the moment
        {
            _phi->FillBoundary(geom->periodicity());
        }

        int nComp() const {return _phi->nComp();}

        bool isReal() const {return _isReal;}
        

        void save(const std::string & folder);


        private:

        multiFab * _phi;
        const geometry * geom;
        bool _isReal;
        std::string _name;
    };



    struct  wavefunctionRegion
    {
        wavefunctionRegion( box boxRegion_,fab * phi ,const geometry* geom_ ) :
        boxRegion(boxRegion_),
        _phi( phi),
        geom(geom_)
         {
             _phiArr=_phi->array();
         }
         
        int minIndex(int d) const {return boxRegion.loVect()[d]; }
        int maxIndex(int d) const {return boxRegion.hiVect()[d];}

       
       






        smallVector minIndex();
        smallVector maxIndex();

        const auto & getBox(){return boxRegion;}
        
        auto & getPhiAsFab(){return *_phi;}
        
        int  nComp() const  {return _phi->nComp(); }

        // enable get access to wavefunction data
        template<class T, std::enable_if_t<std::is_same<T,fab>::value > * = nullptr   >
        auto & getPhi() { return getPhiAsFab();}

        template<class T, std::enable_if_t<std::is_same<T,array_t>::value > * = nullptr  >
        const auto & getPhi() { return _phiArr;}


        const auto & getGeometry() const {return *geom;};
        
        auto size( int d) const
        {
            return maxIndex(d) - minIndex(d) + 1;
        }

        

        const auto size() const
        {
            auto tmpSize = size(0);

            for(int d=1;d<DIMENSIONS;d++)
            {
                tmpSize*=size(d);
            }

            return tmpSize*nComp();
        }

    private:
        box boxRegion;
        fab * _phi;
        const geometry * geom;
        gp::array_t _phiArr;
    };


struct  constWavefunctionRegion
    {
        constWavefunctionRegion( box boxRegion_,const fab * phi ,const geometry* geom_ ) :
        boxRegion(boxRegion_),
        _phi( phi),
        geom(geom_)
         {
             _phiArr=_phi->const_array();
         }

        int minIndex(int d) const {return boxRegion.loVect()[d]; }
        int maxIndex(int d) const {return boxRegion.hiVect()[d];}
        

        auto size( int d) const
        {
            return minIndex(d) - maxIndex(d) + 1;
        }

        const auto size() const
        {
            auto tmpSize = size(0);

            for(int d=1;d<DIMENSIONS;d++)
            {
                tmpSize*=size(d);
            }

            return tmpSize*nComp();
        }




        int  nComp() const  {return _phi->nComp(); }


        smallVector minIndex() const ;
        smallVector maxIndex() const ;

        const auto & getBox(){return boxRegion;}
        
        const auto & getPhiAsFab() const {return *_phi;}


        template<class T, std::enable_if_t<std::is_same<T,fab>::value > * = nullptr   >
        const auto & getPhi() const { return getPhiAsFab();}

        template<class T, std::enable_if_t<std::is_same<T,array_t>::value > * = nullptr  >
        const auto & getPhi() const { return _phiArr;}


        const auto & getGeometry() const {return *geom;};
        
    private:
        box boxRegion;
        const fab * _phi;
        const geometry * geom;
        gp::constArray_t _phiArr;


    };

};

#endif