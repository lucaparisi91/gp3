
#include "AmrLevelDiff.h"
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_StateDescriptor.H>


using namespace amrex;

int      AmrLevelDiff::verbose         = 0;
Real     AmrLevelDiff::cfl             = 0.9;

int      AmrLevelDiff::NUM_STATE       = 1;  // One variable in the state
int      AmrLevelDiff::NUM_GROW        = 3;  // number of ghost cells


void nullfill(Real *data, const int *lo, const int *hi, const int *dom_lo, const int *dom_hi, const Real *dx, const Real *grd_lo, const Real *time, const int *bc)
{

}
AmrLevelDiff::AmrLevelDiff()
{
    
}

//
//The basic constructor.
//
AmrLevelDiff::AmrLevelDiff (Amr&            papa,
     	                  int             lev,
                          const Geometry& level_geom,
                          const BoxArray& bl,
                          const DistributionMapping& dm,
                          Real            time)
    :
    AmrLevel(papa,lev,level_geom,bl,dm,time) 
{
    std::cout << "Built level" << std::endl;
}


AmrLevelDiff::~AmrLevelDiff () 
{
    
}

//
//Restart from a checkpoint file.
//
void
AmrLevelDiff::restart (Amr&          papa,
	              std::istream& is,
                      bool          bReadSpecial)
{
    AmrLevel::restart(papa,is,bReadSpecial);

}

void 
AmrLevelDiff::checkPoint (const std::string& dir,
		         std::ostream&      os,
                         VisMF::How         how,
                         bool               dump_old) 
{
  AmrLevel::checkPoint(dir, os, how, dump_old);
#ifdef AMREX_PARTICLES
  if (do_tracers and level == 0) {
    TracerPC->Checkpoint(dir, "Tracer", true);
  }
#endif
}

//
//Write a plotfile to specified directory.
//
void
AmrLevelDiff::writePlotFile (const std::string& dir,
	 	            std::ostream&      os,
                            VisMF::How         how)
{

    AmrLevel::writePlotFile (dir,os,how);
}

//
//Define data descriptors.
//

void
AmrLevelDiff::variableSetUp ()
{
    BL_ASSERT(desc_lst.size() == 0);

    desc_lst.addDescriptor(Phi_Type,IndexType::TheCellType(),
                           StateDescriptor::Point,0,NUM_STATE,
			   & quadratic_interp);

    int lo_bc[BL_SPACEDIM];
    int hi_bc[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; ++i) {
	lo_bc[i] = hi_bc[i] = BCType::int_dir;   // periodic boundaries
    }
    
    BCRec bc(lo_bc, hi_bc);

    desc_lst.setComponent(Phi_Type, 0, "phi", bc, 
			  StateDescriptor::BndryFunc(nullfill));
}

//
//Cleanup data descriptors at end of run.
//
void
AmrLevelDiff::variableCleanUp () 
{
    desc_lst.clear();
#ifdef AMREX_PARTICLES
    TracerPC.reset();
#endif
}

//
//Initialize grid data at problem start-up.
//
void
AmrLevelDiff::initData ()
{
    //
    // Loop over grids, call FORTRAN function to init with data.
    //
    const Real* dx  = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();
    MultiFab& S_new = get_new_data(Phi_Type);
    Real cur_time   = state[Phi_Type].curTime();

    if (verbose) {
        amrex::Print() << "Initializing the data at level " << level << std::endl;
    }

    Real alpha=1;


    for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
    {
        const Box& box     = mfi.validbox();
        const int* lo      = box.loVect();
        const int* hi      = box.hiVect();


        Array4<  Real> const & dataArray = S_new[mfi].array();


        for(int n=0;n<NUM_STATE;n++)
            for (int k=lo[2];k<=hi[2];k++)
                for (int j=lo[1];j<=hi[1];j++) 
                    for (int i=lo[0];i<=hi[0];i++)
                        {
                            Real x =  prob_lo[0] + (i+0.5)*dx[0];
                            Real y =  prob_lo[1] + (i+0.5)*dx[1];
                            Real z =  prob_lo[2] + (i+0.5)*dx[2];

                            Real r2=x*x + y*y + z*z;
                            dataArray(i,j,k,n)=exp(-alpha*r2);
                        }

    }

    if (verbose) {
	amrex::Print() << "Done initializing the level " << level 
                       << " data " << std::endl;
    }
}

//
//Initialize data on this level from another AmrLevelAdv (during regrid).
//
void
AmrLevelDiff::init (AmrLevel &old)
{
    AmrLevelDiff* oldlev = (AmrLevelDiff*) &old;
    //
    // Create new grid data by fillpatching from old.
    //
    Real dt_new    = parent->dtLevel(level);
    Real cur_time  = oldlev->state[Phi_Type].curTime();
    Real prev_time = oldlev->state[Phi_Type].prevTime();
    Real dt_old    = cur_time - prev_time;
    setTimeLevel(cur_time,dt_old,dt_new);

    MultiFab& S_new = get_new_data(Phi_Type);

    FillPatch(old, S_new, 0, cur_time, Phi_Type, 0, NUM_STATE);
}

//
//Initialize data on this level after regridding if old level did not previously exist
//
void
AmrLevelDiff::init ()
{
    Real dt        = parent->dtLevel(level);
    Real cur_time  = getLevel(level-1).state[Phi_Type].curTime();
    Real prev_time = getLevel(level-1).state[Phi_Type].prevTime();

    Real dt_old = (cur_time - prev_time)/(Real)parent->MaxRefRatio(level-1);

    setTimeLevel(cur_time,dt_old,dt);
    MultiFab& S_new = get_new_data(Phi_Type);
    FillCoarsePatch(S_new, 0, cur_time, Phi_Type, 0, NUM_STATE);
}

//
//Advance grids at this level in time.
//
Real
AmrLevelDiff::advance (Real time,
                      Real dt,
                      int  iteration,
                      int  ncycle)
{
    MultiFab& S_mm = get_new_data(Phi_Type);
    Real maxval = S_mm.max(0);
    Real minval = S_mm.min(0);

    amrex::Print() << "phi max = " << maxval << ", min = " << minval  << std::endl;
    for (int k = 0; k < NUM_STATE_TYPE; k++) {
        state[k].allocOldData();
        state[k].swapTimeLevels(dt);
    }

    MultiFab& S_new = get_new_data(Phi_Type);

    const Real prev_time = state[Phi_Type].prevTime();
    const Real cur_time = state[Phi_Type].curTime();
    const Real ctr_time = 0.5*(prev_time + cur_time);

    const Real* dx = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();



    // State with ghost cells
    MultiFab Sborder(grids, dmap, NUM_STATE, NUM_GROW);
    FillPatch(*this, Sborder, NUM_GROW, time, Phi_Type, 0, NUM_STATE);

    {

     const Real* dx  = geom.CellSize();
     const Real* prob_lo = geom.ProbLo();


	for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
	{
	    const Box& bx = mfi.tilebox();
        const int* lo = bx.loVect();
        const int *hi= bx.hiVect();

	    Array4<  Real> const & statein = Sborder[mfi].array();
	    Array4<  Real> const & stateout = S_new[mfi].array();
        for (int k=lo[2];k<=hi[2];k++)
            for (int j=lo[1];j<=hi[1];j++) 
                for (int i=lo[0];i<=hi[0];i++)
                    {
                        Real x =  prob_lo[0] + (i+0.5)*dx[0];
                        Real y =  prob_lo[1] + (j+0.5)*dx[1];
                        Real z =  prob_lo[2] + (k+0.5)*dx[2];

                        Real r2=x*x + y*y + z*z;
                        stateout(i,j,k)=
                        ( statein(i+1,j,k) - statein(i-1,j,k) )/(2*dx[0]) +
                        ( statein(i,j+1,k) - statein(i,j-1,k) )/(2*dx[1])+
                        ( statein(i,j,k+1) - statein(i,j,k-1) )/(2*dx[2])
                        ;
                    }


	}
    }

    return dt;
}



//
//Compute initial time step.
//
Real
AmrLevelDiff::initialTimeStep ()
{
    return 1e-3;
}

//
//Compute initial `dt'.
//
void
AmrLevelDiff::computeInitialDt (int                   finest_level,
	  	               int                   sub_cycle,
                               Vector<int>&           n_cycle,
                               const Vector<IntVect>& ref_ratio,
                               Vector<Real>&          dt_level,
                               Real                  stop_time)
{
    //
    // Grids have been constructed, compute dt for all levels.
    //
    if (level > 0)
        return;

    Real dt_0 = 1.0e+100;
    int n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        dt_level[i] = getLevel(i).initialTimeStep();
    }

}

//
//Compute new `dt'.
//
void
AmrLevelDiff::computeNewDt (int                   finest_level,
		           int                   sub_cycle,
                           Vector<int>&           n_cycle,
                           const Vector<IntVect>& ref_ratio,
                           Vector<Real>&          dt_min,
                           Vector<Real>&          dt_level,
                           Real                  stop_time,
                           int                   post_regrid_flag)
{
    //
    // Set the same time step at all levels
    if (level > 0)
        return;

    for (int i = 0; i <= finest_level; i++)
    {
        dt_level[i] = getLevel(i).initialTimeStep();
    }

}

//
//Do work after timestep().
//
void
AmrLevelDiff::post_timestep (int iteration)
{


    if (level < parent->finestLevel())
        avgDown();

}

//
//Do work after regrid().
//
void
AmrLevelDiff::post_regrid (int lbase, int new_finest) {


}

//
//Do work after a restart().
//
void
AmrLevelDiff::post_restart() 
{

}

//
//Do work after init().
//
void
AmrLevelDiff::post_init (Real stop_time)
{
    if (level > 0)
        return;
    //
    // Average data down from finer levels
    // so that conserved data is consistent between levels.
    //
    int finest_level = parent->finestLevel();
    for (int k = finest_level-1; k>= 0; k--)
        getLevel(k).avgDown();
}

//
//Error estimation for regridding.
//
void
AmrLevelDiff::errorEst (TagBoxArray& tags,
	               int          clearval,
                       int          tagval,
                       Real         time,
                       int          n_error_buf,
                       int          ngrow)
{
    const Real* dx        = geom.CellSize();
    const Real* prob_lo   = geom.ProbLo();

    MultiFab& S_new = get_new_data(Phi_Type);
    tags.setVal(clearval);


    {
        
        for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
        {
            const Box&  tilebx  = mfi.tilebox();

            Array4<  Real> const & state = S_new[mfi].array();
            //Array4<  int> const & tags = tags[mfi].array();
            

        }
    }
}

void AmrLevelDiff::avgDown ()
{
    if (level == parent->finestLevel()) return;
    avgDown(Phi_Type);
}

void
AmrLevelDiff::avgDown (int state_indx)
{
    if (level == parent->finestLevel()) return;

    AmrLevelDiff& fine_lev = getLevel(level+1);
    MultiFab&  S_fine   = fine_lev.get_new_data(state_indx);
    MultiFab&  S_crse   = get_new_data(state_indx);
    
    amrex::average_down(S_fine,S_crse,
                         fine_lev.geom,geom,
                         0,S_fine.nComp(),parent->refRatio(level));
}
