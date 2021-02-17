

namespace gp
{

    using box=amrex::box;
    using geometry=amrex::geometry;
    
    using multiFab = amrex::Multifab;
    using fab = amrex::FAB;
    using dims = amrex::spaceDim;
    using array_t = amrex::Array4<Real> const;
    using const_array_t = amrex::Array4<const Real>; 
    using smallVector = std::array<int,dims>;
    

    struct wavefunction
    {
        // contains the wavefunction at a single level of refinement

        wavefunction(multifab * newWavefunction, multifab * oldWavefunction.geometry * geom);

        const auto & getGeometry() {return *geom;};

        const auto & getPhi() const {return *phi;}

        auto & getPhi() {return *phi;}

        const auto & getGeometry(){return geom;}


        private:

        multiFab * newWavefunction;
        multiFab * oldWavefunction;
        geometry * geom;

    };

    struct  wavefunctionRegion
    {
        // launch region of the wavefunction
        wavefunctionRegion(box* boxRegion_,fab* currentRegionFab_,geometry* geom);



        
        int minIndex(int d) {return box.loVect()[d]; }
        int maxIndex(int d){return box.hiVect()[d];}

        smallVector minIndex();
        smallVector maxIndex();

        const auto & getBox() const {return *currentRegionBox;}
        const auto & getFab() const {return *currentRegionFab;}


        // enable get access to wavefunction data
        template<class T>
        auto & get() {static_assert(std::is_same<T,fab>::value);return currentRegionFab();}

    private:
        box * boxRegion;
        fab * currentRegionFab;
        geometry * geom;

    };


};

