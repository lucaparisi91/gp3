#include "tools.h"
#include "gpExceptions.h"

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>


Real normCartesian( const MultiFab & phi_real , const MultiFab & phi_imag,  const Geometry & geom, int component)
{
    const Real* dx = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();

    Real norm2=0;
    for ( MFIter mfi( phi_real); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.validbox();
        const int* lo = bx.loVect();
        const int *hi= bx.hiVect();
        Array4< const Real> const & phi_real_box = phi_real[mfi].const_array();
        Array4< const Real> const & phi_imag_box = phi_imag[mfi].const_array();
        const int j=0;
        const int k=0;

#if AMREX_SPACEDIM == 1
            for (int i=lo[0];i<=hi[0];i++)
            {
             
                norm2+=(phi_real_box(i,j,0,component)*phi_real_box(i,j,0,component)
                    + phi_imag_box(i,j,0,component)*phi_imag_box(i,j,0,component));
            }
    
#endif


#if AMREX_SPACEDIM == 2
        for (int j=lo[1];j<=hi[1];j++)
            for (int i=lo[0];i<=hi[0];i++)
            {
             
                norm2+=(phi_real_box(i,j,0,component)*phi_real_box(i,j,0,component)
                    + phi_imag_box(i,j,0,component)*phi_imag_box(i,j,0,component));
            }
    
#endif




#if AMREX_SPACEDIM == 3
        for (int k=lo[2];k<=hi[2];k++)
            for (int j=lo[1];j<=hi[1];j++)
                for (int i=lo[0];i<=hi[0];i++)
                    {

                        norm2+=(phi_real_box(i,j,k,component)*phi_real_box(i,j,k,component)
                        + phi_imag_box(i,j,k,component)*phi_imag_box(i,j,k,component));
                    }
#endif      
    }

    amrex::ParallelAllReduce::Sum(norm2,MPI_COMM_WORLD);

    Real dV=1;

    for (int d=0;d< AMREX_SPACEDIM;d++)
    {
        dV*=dx[d];
    }

    return norm2*dV;
}


Real normSphericalSymmetry( const MultiFab & phi_real , const MultiFab & phi_imag,  const Geometry & geom, int component)
{
    const Real* dx = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();

    Real norm2=0;
    for ( MFIter mfi( phi_real); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.validbox();
        const int* lo = bx.loVect();
        const int *hi= bx.hiVect();
        Array4< const Real> const & phi_real_box = phi_real[mfi].const_array();
        Array4< const Real> const & phi_imag_box = phi_imag[mfi].const_array();
        int j=0;
        int k=0;

        for (int i=lo[0];i<=hi[0];i++)
            {
                Real r =  prob_lo[0] + (i+0.5)*dx[0];

                norm2+=(phi_real_box(i,j,k,component)*phi_real_box(i,j,k,component)
                        + phi_imag_box(i,j,k,component)*phi_imag_box(i,j,k,component)) * r*r;
                    }

    }

    amrex::ParallelAllReduce::Sum(norm2,MPI_COMM_WORLD);


    Real dV=1;

    for (int d=0;d<amrex::SpaceDim;d++)
    {
        dV*=dx[d];
    }

    return norm2*dV*4*M_PI;
}

Real norm( const MultiFab & phi_real , const MultiFab & phi_imag,  const Geometry & geom, int component)
{
    if (geom.Coord() == 0) 
    {
        return normCartesian(phi_real,phi_imag,geom,component);
    }
    else if(geom.Coord()==2)
    {
        return normSphericalSymmetry(phi_real,phi_imag,geom,component);
    }
    else
    {
        throw invalidInput("Unkown coordinates in normalization function.");
        
    }
    
}



void normalize(  MultiFab & phi_real ,  MultiFab & phi_imag,  const Geometry & geom, Real N)
{
    /* Normalizes each component indipendently  */
    for (int i=0;i<phi_real.nComp();i++)
        {
            Real norm2=norm(phi_real,phi_imag,geom,i);
            
            
           //std::cout << "norm:  " << norm2 << std::endl;

            phi_real.mult( std::sqrt(N)/std::sqrt(norm2),i,1);
            phi_imag.mult(std::sqrt(N)/std::sqrt(norm2),i,1);
    }
}



std::tuple< BoxArray , Geometry , DistributionMapping  , std::array<BC,AMREX_SPACEDIM>,   std::array<BC,AMREX_SPACEDIM>  >

createGeometry( const json_t & settings)
{
    BoxArray ba;
    Geometry geom;
    Vector<int> is_periodic(AMREX_SPACEDIM,1);

    std::array<size_t,AMREX_SPACEDIM> shape;
    std::array<Real,AMREX_SPACEDIM> lower_edges;
    std::array<Real,AMREX_SPACEDIM> higher_edges;
    
    std::string coordinates = settings["coordinates"].get<std::string>();


    for (int i=0;i<AMREX_SPACEDIM;i++)
    {
        lower_edges[i]=settings["domain"][i][0].get<Real>();
        higher_edges[i]=settings["domain"][i][1].get<Real>();
        shape[i] = settings["shape"][i].get<int>();

    }

    IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
    IntVect dom_hi(AMREX_D_DECL( shape[0]-1, shape[1]-1, shape[2]-1));
    Box domain(dom_lo, dom_hi);
    ba.define(domain);

    if (settings.contains("maxGridSize") )
    {
        int max_grid_size = settings["maxGridSize"].get<int>();

        ba.maxSize(max_grid_size);

    }
   
    RealBox real_box({AMREX_D_DECL( lower_edges[0],lower_edges[1],lower_edges[2]) },
                         {AMREX_D_DECL( higher_edges[0], higher_edges[1], higher_edges[2] )});
                        
    int coord = 0;

    if (coordinates == "cartesian")
    {
        coord=0;
    }
    else if (coordinates == "spherical")
    {
        coord=2;

        is_periodic[0]=0;

    }
    else
    {
        throw missingImplementation("Unkown coordinates :" + coordinates );
    }

    std::array<BC,AMREX_SPACEDIM> low_bc , high_bc ;

    if (settings.contains("bc") )
    {
        std::tie(low_bc, high_bc) = readBC(settings["bc"]);

        for (int d=0;d<AMREX_SPACEDIM;d++)
        {
            if ( settings["bc"][d] != "periodic" )
            {
                is_periodic[d]=0;
            }    
        }

    }
    else
    {
        auto bc_settings =
        #if AMREX_SPACEDIM == 1     
        R"( ["periodic"   ])"_json
        #endif
         #if AMREX_SPACEDIM == 3     
        R"( ["periodic" , "periodic" , "periodic"  ])"_json
        #endif
        ;
        std::tie(low_bc, high_bc) = readBC(bc_settings);

    }
    

    geom.define(domain,&real_box,coord,is_periodic.data());

    DistributionMapping dm(ba);

    return {ba, geom, dm, low_bc, high_bc} ;


}

void fill(MultiFab & realState, MultiFab & imagState, py::array_t<std::complex<Real> > initialCondition , Geometry & geom)
{
    auto psi = initialCondition.unchecked<AMREX_SPACEDIM>();

    LOOP(realState,geom)

#if AMREX_SPACEDIM == 3
       data(i,j,k)=std::real( psi(i,j,k) );
#endif

#if AMREX_SPACEDIM == 1
       data(i,j,k)=std::real( psi(i) );

#endif



    ENDLOOP

    LOOP(imagState,geom)

#if AMREX_SPACEDIM == 3
       data(i,j,k)=std::imag( psi(i,j,k) );
#endif

#if AMREX_SPACEDIM == 1
       data(i,j,k)=std::imag( psi(i) );
#endif

    ENDLOOP
}

LinOpBCType toLinOpBCType( const BC & bc )
{

    if ( bc == BC::DRICHLET) 
    {
        return LinOpBCType::Dirichlet;
    }
    else if (bc == BC::NEUMANNN)
    {
        return LinOpBCType::Neumann;
    }
    else
    {
        return LinOpBCType::Periodic;
    }

}


std::tuple<std::array<BC,AMREX_SPACEDIM>, std::array<BC,AMREX_SPACEDIM> > readBC( const json_t & j)
{
    std::array<BC,AMREX_SPACEDIM> low_bc;
    std::array<BC,AMREX_SPACEDIM> high_bc;

    for (int d=0;d<AMREX_SPACEDIM;d++)
    {
        BC bc = BC::PERIODIC;
        auto bc_type = j[d].get<std::string>(); 
        if ( bc_type == "neumann" )
        {
            bc=BC::NEUMANNN;
        }
        else
        { if (bc_type == "drichlet")
            {
                bc=BC::DRICHLET;
            }
            else if (bc_type != "periodic")
            {
                throw invalidInput("Unkown bc " + bc_type);
            }
        }
        low_bc[d]=bc;
        high_bc[d] = bc;
    }

    return std::tuple(low_bc,high_bc);
}


auto toMultiFabBC(const BC & bc)
{
    if (bc== BC::PERIODIC)
    {
        return BCType::int_dir;
    }
    else if (bc == BC::NEUMANNN)
    {
        return BCType::reflect_even;
    }
    else if (bc == BC::DRICHLET)
    {
        return BCType::ext_dir;
    }
    else
    {
        throw invalidInput("Unkown boundary condition");
    }

}


BCRec toMultiFabBC(const std::array<BC,AMREX_SPACEDIM> & bc_low, const std::array<BC,AMREX_SPACEDIM> & bc_high  )
{
  BCRec bc;
    
    for (int d=0;d<AMREX_SPACEDIM;d++)
    {  
        bc.setLo(d, toMultiFabBC(bc_low[d]) ); 
        bc.setHi(d, toMultiFabBC(bc_high[d]));
    }

    return bc;
}


