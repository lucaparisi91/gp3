
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_MLPoisson.H>
#include <AMReX_MLMG.H>

#ifdef AMREX_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

#include "AmrCoreDiff.h"

using namespace amrex;

// constructor - reads in parameters from inputs file
//             - sizes multilevel arrays and data structures
//             - initializes BCRec boundary condition object
AmrCoreDiff::AmrCoreDiff(const RealBox* rb, int max_level_in,
             const Vector<int>& n_cell_in, int coord) : AmrCore::AmrCore(rb,max_level_in,n_cell_in,coord)
{
    //ReadParameters();

    // Geometry on all levels has been defined already.

    // No valid BoxArray and DistributionMapping have been defined.
    // But the arrays for them have been resized.

    int nlevs_max = max_level + 1;


    istep.resize(nlevs_max, 0);
    nsubsteps.resize(nlevs_max, 1);
    if (do_subcycle) {
        for (int lev = 1; lev <= max_level; ++lev) {
            nsubsteps[lev] = MaxRefRatio(lev-1);
        }
    }

    t_new.resize(nlevs_max, 0.0);
    t_old.resize(nlevs_max, -1.e100);
    dt.resize(nlevs_max, 1.e100);

    phi_new.resize(nlevs_max);
    phi_old.resize(nlevs_max);

    flux_new.resize(nlevs_max);
    flux_old.resize(nlevs_max);




    facevel.resize(nlevs_max);

    // periodic boundaries
    int bc_lo[] = {BCType::int_dir, BCType::int_dir, BCType::int_dir};
    int bc_hi[] = {BCType::int_dir, BCType::int_dir, BCType::int_dir};

/*
    // walls (Neumann)
    int bc_lo[] = {FOEXTRAP, FOEXTRAP, FOEXTRAP};
    int bc_hi[] = {FOEXTRAP, FOEXTRAP, FOEXTRAP};
*/

    bcs.resize(1);     // Setup 1-component
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        // lo-side BCs
        if (bc_lo[idim] == BCType::int_dir  ||  // periodic uses "internal Dirichlet"
            bc_lo[idim] == BCType::foextrap ||  // first-order extrapolation
            bc_lo[idim] == BCType::ext_dir ) {  // external Dirichlet
            bcs[0].setLo(idim, bc_lo[idim]);
        }
        else {
            amrex::Abort("Invalid bc_lo");
        }

        // hi-side BCSs
        if (bc_hi[idim] == BCType::int_dir  ||  // periodic uses "internal Dirichlet"
            bc_hi[idim] == BCType::foextrap ||  // first-order extrapolation
            bc_hi[idim] == BCType::ext_dir ) {  // external Dirichlet
            bcs[0].setHi(idim, bc_hi[idim]);
        }
        else {
            amrex::Abort("Invalid bc_hi");
        }
    }

    // stores fluxes at coarse-fine interface for synchronization
    // this will be sized "nlevs_max+1"
    // NOTE: the flux register associated with flux_reg[lev] is associated
    // with the lev/lev-1 interface (and has grid spacing associated with lev-1)
    // therefore flux_reg[0] is never actually used in the reflux operation
    flux_reg.resize(nlevs_max+1);
    regrid_int=10 ;
    plot_int = 100;
    max_step=10000;

    blocking_factor={ };
    n_error_buf={};

    for (int i=0;i<nlevs_max + 1 ;i++)
    {
        blocking_factor.push_back({1,1,1});
        n_error_buf.push_back({0,0,0});

    }

    for(int i=0;i<nlevs_max;i++)
    {
        geom[i].setPeriodicity({1,1,1});
    }

    advanceMethod = advanceMethod_t::FillPatchAdvance;

}

AmrCoreDiff::~AmrCoreDiff ()
{
}

// advance solution to final time
void
AmrCoreDiff::Evolve ()
{
    Real cur_time = t_new[0];
    int last_plot_file_step = 0;

    for (int step = istep[0]; step < max_step && cur_time < stop_time; ++step)
    {
        amrex::Print() << "\nCoarse STEP " << step+1 << " starts ..." << std::endl;

        ComputeDt();

        int lev = 0;
        int iteration = 1;
        if (do_subcycle)
            timeStepWithSubcycling(lev, cur_time, iteration);
        else
            timeStepNoSubcycling(cur_time, iteration);

        cur_time += dt[0];

        // sum phi to check conservation
        Real sum_phi = phi_new[0].sum();

        amrex::Print() << "Coarse STEP " << step+1 << " ends." << " TIME = " << cur_time
                       << " DT = " << dt[0] << " Sum(Phi) = " << sum_phi << std::endl;

        // sync up time
        for (lev = 0; lev <= finest_level; ++lev) {
            t_new[lev] = cur_time;
        }

        if (plot_int > 0 && (step+1) % plot_int == 0) {
            last_plot_file_step = step+1;
            WritePlotFile();
        }

        if (chk_int > 0 && (step+1) % chk_int == 0) {
            WriteCheckpointFile();
        }

#ifdef AMREX_MEM_PROFILING
        {
            std::ostringstream ss;
            ss << "[STEP " << step+1 << "]";
            MemProfiler::report(ss.str());
        }
#endif

        if (cur_time >= stop_time - 1.e-6*dt[0]) break;
    }

    if (plot_int > 0 && istep[0] > last_plot_file_step) {
        WritePlotFile();
    }
}

// initializes multilevel data
void
AmrCoreDiff::InitData ()
{
    if (restart_chkfile == "") {
        // start simulation from the beginning
        const Real time = 0.0;
        InitFromScratch(time);
        AverageDown();

        if (chk_int > 0) {
            WriteCheckpointFile();
        }

    }
    else {
        // restart from a checkpoint
        ReadCheckpointFile();
    }

    if (plot_int > 0) {
        WritePlotFile();
    }
}

// Make a new level using provided BoxArray and DistributionMapping and 
// fill with interpolated coarse level data.
// overrides the pure virtual function in AmrCore
void
AmrCoreDiff::MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba,
				    const DistributionMapping& dm)
{
    const int ncomp = phi_new[lev-1].nComp();
    const int nghost = phi_new[lev-1].nGrow();
    
    phi_new[lev].define(ba, dm, ncomp, nghost);
    phi_old[lev].define(ba, dm, ncomp, nghost);

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    // This clears the old MultiFab and allocates the new one
    for (int idim = 0; idim < AMREX_SPACEDIM; idim++)
    {
	facevel[lev][idim] = MultiFab(amrex::convert(ba,IntVect::TheDimensionVector(idim)), dm, 1, 1);
    }

    if (lev > 0 && do_reflux) {
	flux_reg[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, ncomp));
    }

    FillCoarsePatch(lev, time, phi_new[lev], 0, ncomp);
}

// Remake an existing level using provided BoxArray and DistributionMapping and 
// fill with existing fine and coarse data.
// overrides the pure virtual function in AmrCore
void
AmrCoreDiff::RemakeLevel (int lev, Real time, const BoxArray& ba,
			 const DistributionMapping& dm)
{
    const int ncomp = phi_new[lev].nComp();
    const int nghost = phi_new[lev].nGrow();

    MultiFab new_state(ba, dm, ncomp, nghost);
    MultiFab old_state(ba, dm, ncomp, nghost);

    FillPatch(lev, time, new_state, 0, ncomp);

    std::swap(new_state, phi_new[lev]);
    std::swap(old_state, phi_old[lev]);

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    // This clears the old MultiFab and allocates the new one
    for (int idim = 0; idim < AMREX_SPACEDIM; idim++)
    {
	facevel[lev][idim] = MultiFab(amrex::convert(ba,IntVect::TheDimensionVector(idim)), dm, 1, 1);
    }

    if (lev > 0 && do_reflux) {
	flux_reg[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, ncomp));
    }    
}

// Delete level data
// overrides the pure virtual function in AmrCore
void
AmrCoreDiff::ClearLevel (int lev)
{
    phi_new[lev].clear();
    phi_old[lev].clear();
    flux_reg[lev].reset(nullptr);
}

// Make a new level from scratch using provided BoxArray and DistributionMapping.
// Only used during initialization.
// overrides the pure virtual function in AmrCore
void AmrCoreDiff::MakeNewLevelFromScratch (int lev, Real time, const BoxArray& ba,
					  const DistributionMapping& dm)
{
    const int ncomp = 1;
    const int nghost = 3;


    phi_new[lev].define(ba, dm, ncomp, nghost);
    phi_old[lev].define(ba, dm, ncomp, nghost);


    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    // This clears the old MultiFab and allocates the new one
    for (int idim = 0; idim < AMREX_SPACEDIM; idim++)
    {
	facevel[lev][idim] = MultiFab(amrex::convert(ba,IntVect::TheDimensionVector(idim)), dm, 1, 1);
    flux_new[lev][idim] = MultiFab(amrex::convert(ba,IntVect::TheDimensionVector(idim)), dm, 1, 0);
    flux_old[lev][idim] = MultiFab(amrex::convert(ba,IntVect::TheDimensionVector(idim)), dm, 1, 0);
    
    }

    if (lev > 0 && do_reflux) {
	flux_reg[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, ncomp));
    }

    MultiFab& state = phi_new[lev];

    const auto problo = Geom(lev).ProbLoArray();
    const auto dx     = Geom(lev).CellSizeArray();
    int nComps = state.nComp();


    Real var=1;
    Real alpha=1/(2.*var);


#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(state,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Array4<Real> fab = state[mfi].array();
        const Box& box = mfi.tilebox();
        const int* lo      = box.loVect();
        const int* hi      = box.hiVect();

        for(int n=0;n<nComps;n++)
                for (int k=lo[2];k<=hi[2];k++)
                    for (int j=lo[1];j<=hi[1];j++) 
                        for (int i=lo[0];i<=hi[0];i++)
                            {
                                Real x =  problo[0] + (i+0.5)*dx[0];
                                Real y =  problo[1] + (j+0.5)*dx[1];
                                Real z =  problo[2] + (k+0.5)*dx[2];

                                Real r2=x*x + y*y + z*z;
                                
                                    fab(i,j,k,n)=exp(-alpha*r2);
                                
                                
                            }

    }
}

// tag all cells for refinement
// overrides the pure virtual function in AmrCore
void
AmrCoreDiff::ErrorEst (int lev, TagBoxArray& tags, Real /*time*/, int /*ngrow*/)
{


//    const int clearval = TagBox::CLEAR;
    const int   tagval = TagBox::SET;

    const MultiFab& state = phi_new[lev];
    const auto problo = Geom(lev).ProbLoArray();
    const auto dx     = Geom(lev).CellSizeArray();

    
    std::vector<Real> phi_lims{0.01,0.1,0.5,0.9} ;

#ifdef _OPENMP
#pragma omp parallel if(Gpu::notInLaunchRegion())
#endif
    {
	
	for (MFIter mfi(state); mfi.isValid(); ++mfi)
	{
	    const Box& bx  = mfi.validbox();
        const auto statefab = state.array(mfi);
        const auto tagfab  = tags.array(mfi);
        const int* lo      = bx.loVect();
        const int* hi      = bx.hiVect();


        for(int n=0;n<state.nComp();n++)
          for (int k=lo[2];k<=hi[2];k++)
                for (int j=lo[1];j<=hi[1];j++) 
                    for (int i=lo[0];i<=hi[0];i++)
                        {

                            Real x =  problo[0] + (i+0.5)*dx[0];
                            Real y =  problo[1] + (j+0.5)*dx[1];
                            Real z =  problo[2] + (k+0.5)*dx[2];

                            Real r2=x*x + y*y + z*z;
/*                              if ( (i>=10 and i<=20) and
                                (j>=10 and j<=20) and 
                                (k>=10 and k<=20) )  */
                            if (statefab(i,j,k,n)>=phi_lims[lev])    
                            {
                                tagfab(i,j,k,n)=tagval;
                            }
                            else
                            {
                                tagfab(i,j,k,n)=TagBox::CLEAR;
                            }
                            
                        }

    }
}
}

// read in some parameters from inputs file
void
AmrCoreDiff::ReadParameters ()
{
    {
	ParmParse pp;  // Traditionally, max_step and stop_time do not have prefix.
	pp.query("max_step", max_step);
	pp.query("stop_time", stop_time);
    }

    {
	ParmParse pp("amr"); // Traditionally, these have prefix, amr.

	pp.query("regrid_int", regrid_int);
	pp.query("plot_file", plot_file);
	pp.query("plot_int", plot_int);
	pp.query("chk_file", chk_file);
	pp.query("chk_int", chk_int);
        pp.query("restart",restart_chkfile);
    }

    {
	ParmParse pp("adv");
	
	pp.query("cfl", cfl);
        pp.query("do_reflux", do_reflux);
        pp.query("do_subcycle", do_subcycle);
    }
}




// set covered coarse cells to be the average of overlying fine cells
void
AmrCoreDiff::AverageDown ()
{
    for (int lev = finest_level-1; lev >= 0; --lev)
    {
	amrex::average_down(phi_new[lev+1], phi_new[lev],
                            geom[lev+1], geom[lev],
                            0, phi_new[lev].nComp(), refRatio(lev));
    }
}

// more flexible version of AverageDown() that lets you average down across multiple levels
void
AmrCoreDiff::AverageDownTo (int crse_lev)
{
    amrex::average_down(phi_new[crse_lev+1], phi_new[crse_lev],
                        geom[crse_lev+1], geom[crse_lev],
                        0, phi_new[crse_lev].nComp(), refRatio(crse_lev));
}

// compute a new multifab by coping in phi from valid region and filling ghost cells
// works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
void
AmrCoreDiff::FillPatch (int lev, Real time, MultiFab& mf, int icomp, int ncomp)
{
    if (lev == 0)
    {
	Vector<MultiFab*> smf;
	Vector<Real> stime;
	GetData(0, time, smf, stime);

        
        {
            CpuBndryFuncFab bndry_func(nullptr);  // Without EXT_DIR, we can pass a nullptr.
            PhysBCFunct<CpuBndryFuncFab> physbc(geom[lev],bcs,bndry_func);
            amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, icomp, ncomp, 
                                        geom[lev], physbc, 0);
        }
    }
    else
    {
	Vector<MultiFab*> cmf, fmf;
	Vector<Real> ctime, ftime;
	GetData(lev-1, time, cmf, ctime);
	GetData(lev  , time, fmf, ftime);

	Interpolater* mapper = 
    //& cell_cons_interp;
    &quartic_interp;
        
        {
            CpuBndryFuncFab bndry_func(nullptr);  // Without EXT_DIR, we can pass a nullptr.
            PhysBCFunct<CpuBndryFuncFab> cphysbc(geom[lev-1],bcs,bndry_func);
            PhysBCFunct<CpuBndryFuncFab> fphysbc(geom[lev],bcs,bndry_func);

            amrex::FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime,
                                      0, icomp, ncomp, geom[lev-1], geom[lev],
                                      cphysbc, 0, fphysbc, 0, refRatio(lev-1),
                                      mapper, bcs, 0);
        }
    }
}

// fill an entire multifab by interpolating from the coarser level
// this comes into play when a new level of refinement appears
void
AmrCoreDiff::FillCoarsePatch (int lev, Real time, MultiFab& mf, int icomp, int ncomp)
{
    BL_ASSERT(lev > 0);

    Vector<MultiFab*> cmf;
    Vector<Real> ctime;
    GetData(lev-1, time, cmf, ctime);
    Interpolater* mapper = 
    //&cell_cons_interp; 
    &quartic_interp;

    if (cmf.size() != 1) {
	amrex::Abort("FillCoarsePatch: how did this happen?");
    }

    {
        CpuBndryFuncFab bndry_func(nullptr);  // Without EXT_DIR, we can pass a nullptr.
        PhysBCFunct<CpuBndryFuncFab> cphysbc(geom[lev-1],bcs,bndry_func);
        PhysBCFunct<CpuBndryFuncFab> fphysbc(geom[lev],bcs,bndry_func);

        amrex::InterpFromCoarseLevel(mf, time, *cmf[0], 0, icomp, ncomp, geom[lev-1], geom[lev],
                                     cphysbc, 0, fphysbc, 0, refRatio(lev-1),
                                     mapper, bcs, 0);
    }
}

// utility to copy in data from phi_old and/or phi_new into another multifab
void
AmrCoreDiff::GetData (int lev, Real time, Vector<MultiFab*>& data, Vector<Real>& datatime)
{
    data.clear();
    datatime.clear();

    const Real teps = (t_new[lev] - t_old[lev]) * 1.e-3;

    if (time > t_new[lev] - teps && time < t_new[lev] + teps)
    {
	data.push_back(&phi_new[lev]);
	datatime.push_back(t_new[lev]);
    }
    else if (time > t_old[lev] - teps && time < t_old[lev] + teps)
    {
	data.push_back(&phi_old[lev]);
	datatime.push_back(t_old[lev]);
    }
    else
    {
	data.push_back(&phi_old[lev]);
	data.push_back(&phi_new[lev]);
	datatime.push_back(t_old[lev]);
	datatime.push_back(t_new[lev]);
    }
}

void
AmrCoreDiff::GetFlux(int lev, Real time, Vector<MultiFab*>& data, Vector<Real>& datatime, int dir)
{
    data.clear();
    datatime.clear();

    const Real teps = (t_new[lev] - t_old[lev]) * 1.e-3;

    if (time > t_new[lev] - teps && time < t_new[lev] + teps)
    {
	data.push_back(&(flux_new[lev][dir]));
	datatime.push_back(t_new[lev]);
    }
    else if (time > t_old[lev] - teps && time < t_old[lev] + teps)
    {

	data.push_back(&flux_old[lev][dir]);
	datatime.push_back(t_old[lev]);
    }
    else
    {
	data.push_back(&flux_old[lev][dir]);
	data.push_back(&flux_new[lev][dir]);
	datatime.push_back(t_old[lev]);
	datatime.push_back(t_new[lev]);
    }
}


// Advance a level by dt
// (includes a recursive call for finer levels)
void
AmrCoreDiff::timeStepWithSubcycling(int lev, Real time, int iteration)
{

    if (regrid_int > 0)  // We may need to regrid
    {

        // help keep track of whether a level was already regridded
        // from a coarser level call to regrid
        static Vector<int> last_regrid_step(max_level+1, 0);

        // regrid changes level "lev+1" so we don't regrid on max_level
        // also make sure we don't regrid fine levels again if 
        // it was taken care of during a coarser regrid
        if (lev < max_level && istep[lev] > last_regrid_step[lev]) 
        {
            if (istep[lev] % regrid_int == 0)
            {
                // regrid could add newly refine levels (if finest_level < max_level)
                // so we save the previous finest level index
                int old_finest = finest_level; 
                regrid(lev, time);

                // mark that we have regridded this level already
                for (int k = lev; k <= finest_level; ++k) {
                    last_regrid_step[k] = istep[k];
                }

                // if there are newly created levels, set the time step
                for (int k = old_finest+1; k <= finest_level; ++k) {
                    dt[k] = 
                    //dt[k-1] / MaxRefRatio(k-1);
                    dt[k-1];
                     
                }
            }
        }
    }

    if (Verbose()) {
        amrex::Print() << "[Level " << lev << " step " << istep[lev]+1 << "] ";
        amrex::Print() << "ADVANCE with time = " << t_new[lev] 
                       << " dt = " << dt[lev] << std::endl;
    }

    // Advance a single level for a single time step, and update flux registers

    //t_old[lev] = t_new[lev];
    //t_new[lev] += dt[lev];

    Real t_nph = t_old[lev] + 0.5*dt[lev]; 

    
    AdvancePhiAtLevel(lev, time, dt[lev], iteration, nsubsteps[lev]);

    ++istep[lev];

    if (Verbose())
    {
        amrex::Print() << "[Level " << lev << " step " << istep[lev] << "] ";
        amrex::Print() << "Advanced " << CountCells(lev) << " cells" << std::endl;
    }

    if (lev < finest_level)
    {
        // recursive call for next-finer level
        for (int i = 1; i <= nsubsteps[lev+1]; ++i)
        {
            timeStepWithSubcycling(lev+1, time+(i-1)*dt[lev+1], i);
        }

        

        AverageDownTo(lev); // average lev+1 down to lev
    }
    
}

// Advance all the levels with the same dt
void
AmrCoreDiff::timeStepNoSubcycling (Real time, int iteration)
{

}

// a wrapper for EstTimeStep
void
AmrCoreDiff::ComputeDt ()
{
    Vector<Real> dt_tmp(finest_level+1);

    dt[0]=timeStep;
    for (int lev =  1; lev <= finest_level; ++lev)
    {
        //dt[lev] = dt[lev-1]/MaxRefRatio(lev);
        dt[lev]=timeStep;
        nsubsteps[lev]=1;
    }
    
}

// compute dt from CFL considerations
Real
AmrCoreDiff::EstTimeStep (int lev, Real time)
{
    return timeStep;
}

// get plotfile name
std::string
AmrCoreDiff::PlotFileName (int lev) const
{
    return amrex::Concatenate(plot_file, lev, 5);
}

// put together an array of multifabs for writing
Vector<const MultiFab*>
AmrCoreDiff::PlotFileMF () const
{
    Vector<const MultiFab*> r;
    for (int i = 0; i <= finest_level; ++i) {
	r.push_back(&phi_new[i]);
    }
    return r;
}

// set plotfile variable names
Vector<std::string>
AmrCoreDiff::PlotFileVarNames () const
{
    return {"phi"};
}

// write plotfile to disk
void
AmrCoreDiff::WritePlotFile () const
{
    const std::string& plotfilename = PlotFileName(istep[0]);
    const auto& mf = PlotFileMF();
    const auto& varnames = PlotFileVarNames();
    
    amrex::Print() << "Writing plotfile " << plotfilename << "\n";

    amrex::WriteMultiLevelPlotfile(plotfilename, finest_level+1, mf, varnames,
				   Geom(), t_new[0], istep, refRatio());
}

void
AmrCoreDiff::WriteCheckpointFile () const
{

    // chk00010            write a checkpoint file with this root directory
    // chk00010/Header     this contains information you need to save (e.g., finest_level, t_new, etc.) and also
    //                     the BoxArrays at each level
    // chk00010/Level_0/
    // chk00010/Level_1/
    // etc.                these subdirectories will hold the MultiFab data at each level of refinement

    // checkpoint file name, e.g., chk00010
    const std::string& checkpointname = amrex::Concatenate(chk_file,istep[0]);

    amrex::Print() << "Writing checkpoint " << checkpointname << "\n";

    const int nlevels = finest_level+1;

    // ---- prebuild a hierarchy of directories
    // ---- dirName is built first.  if dirName exists, it is renamed.  then build
    // ---- dirName/subDirPrefix_0 .. dirName/subDirPrefix_nlevels-1
    // ---- if callBarrier is true, call ParallelDescriptor::Barrier()
    // ---- after all directories are built
    // ---- ParallelDescriptor::IOProcessor() creates the directories
    amrex::PreBuildDirectorHierarchy(checkpointname, "Level_", nlevels, true);

    // write Header file
   if (ParallelDescriptor::IOProcessor()) {

       std::string HeaderFileName(checkpointname + "/Header");
       VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
       std::ofstream HeaderFile;
       HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
       HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
		                               std::ofstream::trunc |
                                               std::ofstream::binary);
       if( ! HeaderFile.good()) {
           amrex::FileOpenFailed(HeaderFileName);
       }

       HeaderFile.precision(17);

       // write out title line
       HeaderFile << "Checkpoint file for AmrCoreDiff\n";

       // write out finest_level
       HeaderFile << finest_level << "\n";

       // write out array of istep
       for (int i = 0; i < istep.size(); ++i) {
           HeaderFile << istep[i] << " ";
       }
       HeaderFile << "\n";

       // write out array of dt
       for (int i = 0; i < dt.size(); ++i) {
           HeaderFile << dt[i] << " ";
       }
       HeaderFile << "\n";

       // write out array of t_new
       for (int i = 0; i < t_new.size(); ++i) {
           HeaderFile << t_new[i] << " ";
       }
       HeaderFile << "\n";

       // write the BoxArray at each level
       for (int lev = 0; lev <= finest_level; ++lev) {
           boxArray(lev).writeOn(HeaderFile);
           HeaderFile << '\n';
       }
   }

   // write the MultiFab data to, e.g., chk00010/Level_0/
   for (int lev = 0; lev <= finest_level; ++lev) {
       VisMF::Write(phi_new[lev],
                    amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "phi"));
   }

}

namespace {
// utility to skip to next line in Header
void GotoNextLine (std::istream& is)
{
    constexpr std::streamsize bl_ignore_max { 100000 };
    is.ignore(bl_ignore_max, '\n');
}
}

void
AmrCoreDiff::ReadCheckpointFile ()
{

    amrex::Print() << "Restart from checkpoint " << restart_chkfile << "\n";

    // Header
    std::string File(restart_chkfile + "/Header");

    VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);

    std::string line, word;

    // read in title line
    std::getline(is, line);

    // read in finest_level
    is >> finest_level;
    GotoNextLine(is);

    // read in array of istep
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            istep[i++] = std::stoi(word);
        }
    }

    // read in array of dt
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            dt[i++] = std::stod(word);
        }
    }

    // read in array of t_new
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            t_new[i++] = std::stod(word);
        }
    }

    for (int lev = 0; lev <= finest_level; ++lev) {

        // read in level 'lev' BoxArray from Header
        BoxArray ba;
        ba.readFrom(is);
        GotoNextLine(is);

        // create a distribution mapping
        DistributionMapping dm { ba, ParallelDescriptor::NProcs() };

        // set BoxArray grids and DistributionMapping dmap in AMReX_AmrMesh.H class
        SetBoxArray(lev, ba);
        SetDistributionMap(lev, dm);

        // build MultiFab and FluxRegister data
        int ncomp = 1;
        int nghost = 0;
        phi_old[lev].define(grids[lev], dmap[lev], ncomp, nghost);
        phi_new[lev].define(grids[lev], dmap[lev], ncomp, nghost);

        if (lev > 0 && do_reflux) {
            flux_reg[lev].reset(new FluxRegister(grids[lev], dmap[lev], refRatio(lev-1), lev, ncomp));
        }

        // build face velocity MultiFabs
        for (int idim = 0; idim < AMREX_SPACEDIM; idim++)
        {
	    facevel[lev][idim] = MultiFab(amrex::convert(ba,IntVect::TheDimensionVector(idim)), dm, 1, 1);
        }
    }

    // read in the MultiFab data
    for (int lev = 0; lev <= finest_level; ++lev) {
        VisMF::Read(phi_new[lev],
                    amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "phi"));
    }

}


void
AmrCoreDiff::AdvancePhiAtLevel (int lev, Real time, Real dt_lev, int iteration, int ncycle)
{
    if (advanceMethod == advanceMethod_t::FillPatchAdvance )
    {
        AdvancePhiAtLevelFillPatch(lev, time, dt_lev, iteration, ncycle);
    }
    else if        (advanceMethod == advanceMethod_t::LinOpAdvance )
    {
        AdvancePhiAtLevelMLMG(lev, time, dt_lev, iteration, ncycle);
    }   

}
void
AmrCoreDiff::AdvancePhiAtLevelMLMG (int lev, Real time, Real dt_lev, int /*iteration*/, int /*ncycle*/)
{
    constexpr int num_grow = 3;

    std::swap(phi_old[lev], phi_new[lev]);
    t_old[lev] = t_new[lev];
    t_new[lev] += dt_lev;

    const auto dx = geom[lev].CellSizeArray();

    LPInfo info;
    MLPoisson mlpoisson({Geom(lev)}, {grids[lev]}, {dmap[lev]}, info);

    mlpoisson.setDomainBC({AMREX_D_DECL(LinOpBCType::Periodic,
                                                LinOpBCType::Periodic,
                                                LinOpBCType::Periodic)},
                                  {AMREX_D_DECL(LinOpBCType::Periodic,
                                                LinOpBCType::Periodic,
                                                LinOpBCType::Periodic)});

    if (lev>0)
    {
        mlpoisson.setCoarseFineBC(&phi_old[lev-1],MaxRefRatio(lev-1));
    }

    mlpoisson.setLevelBC(0,&phi_old[lev]);

    MLMG mlmg(mlpoisson);

    mlmg.apply({&phi_new[lev]},{&phi_old[lev]});


    //phi_new[lev].mult(dt_lev,num_grow);

    //phi_new[lev].plus(phi_old[lev],0,phi_new[lev].nComp(),num_grow);



    
  /*   // State with ghost cells
    MultiFab Sborder(grids[lev], dmap[lev], phi_new[lev].nComp(), num_grow);
    // State with ghost cells
    MultiFab S(grids[lev], dmap[lev], phi_new[lev].nComp(), 0);
    

     if (lev==1)
    {
        FillCoarsePatch(lev, time, Sborder, 0, phi_new[lev].nComp());
        S.copy(phi_old[lev]);
        Sborder.copy(S);



    }
    else
    {
        FillPatch(lev, time, Sborder, 0, Sborder.nComp());
    
    }
    
    
    

    for (MFIter mfi(Sborder,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Array4<Real> const & fab = Sborder[mfi].array();
        Array4<Real> const & phi_new_fab = phi_new[lev][mfi].array();    
        const Box& box = mfi.tilebox();
        const int* lo      = box.loVect();
        const int* hi      = box.hiVect();
        const auto problo = Geom(lev).ProbLoArray();
        const auto dx     = Geom(lev).CellSizeArray();

        for(int n=0;n<Sborder.nComp();n++)
            for (int k=lo[2];k<=hi[2];k++)
                for (int j=lo[1];j<=hi[1];j++) 
                    for (int i=lo[0];i<=hi[0];i++)
                        {
                            Real x =  problo[0] + (i+0.5)*dx[0];
                            Real y =  problo[1] + (j+0.5)*dx[1];
                            Real z =  problo[2] + (k+0.5)*dx[2];

                            Real r2=x*x + y*y + z*z;

                            std::cout << fab(i,j,k,n) << std::endl;
                            phi_new_fab(i,j,k,n)=
                            //fab(i,j,k,n);
        
        
                            //   std::pow( ( fab(i+1,j,k,n) - fab(i-1,j,k,n))/(dx[0]*2) ,2) +
                            // std::pow( ( fab(i,j+1,k,n) - fab(i,j-1,k,n))/(dx[1] *2),2) +
                            // std::pow( ( fab(i,j,k+1,n) - fab(i,j,k-1,n))/(dx[2]*2) ,2) ; 
                             

                             
                            //   ( fab(i+1,j,k,n) + fab(i-1,j,k,n) -2* fab(i,j,k,n) )/(dx[0]*dx[0] )+
                            // ( fab(i,j+1,k,n) + fab(i,j-1,k,n) -2* fab(i,j,k,n) )/(dx[1]*dx[1]) +
                            // ( fab(i,j,k+1,n) + fab(i,j,k-1,n) -2* fab(i,j,k,n) )/(dx[2]*dx[2]);  

                            ( -1./12. * fab(i-2,j,k,n) + 4./3 *fab(i-1,j,k,n) -5./2* fab(i,j,k,n) + 4/3.*fab(i+1,j,k,n)  -1./12*fab(i+2,j,k,n)  )/(dx[0]*dx[0] ) + 
                            ( -1./12. * fab(i,j-2,k,n) + 4./3 *fab(i,j-1,k,n) -5./2* fab(i,j,k,n) + 4/3.*fab(i,j+1,k,n)  -1./12*fab(i,j+2,k,n)  )/(dx[1]*dx[1] ) +
                            ( -1./12. * fab(i,j,k-2,n) + 4./3 *fab(i,j,k-1,n) -5./2* fab(i,j,k,n) + 4/3.*fab(i,j,k+1,n)  -1./12*fab(i,j,k+2,n)  )/(dx[2]*dx[2] );


                            





                        }
    } */
    
}

// advance all levels for a single time step
void
AmrCoreDiff::AdvancePhiAtLevelFillPatch (int lev, Real time, Real dt_lev, int /*iteration*/, int /*ncycle*/)
{

    constexpr int num_grow = 3;

    std::swap(phi_old[lev], phi_new[lev]);
    t_old[lev] = t_new[lev];
    t_new[lev] += dt_lev;

    const auto dx = geom[lev].CellSizeArray();

    


    
   // State with ghost cells
    MultiFab Sborder(grids[lev], dmap[lev], phi_new[lev].nComp(), num_grow);
    // State with ghost cells
    MultiFab S(grids[lev], dmap[lev], phi_new[lev].nComp(), 0);
    

     if (lev==1)
    {
        FillCoarsePatch(lev, time, Sborder, 0, phi_new[lev].nComp());
        S.copy(phi_old[lev]);
        Sborder.copy(S);



    }
    else
    {
        FillPatch(lev, time, Sborder, 0, Sborder.nComp());
    
    }
    
    
    

    for (MFIter mfi(Sborder,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Array4<Real> const & fab = Sborder[mfi].array();
        Array4<Real> const & phi_new_fab = phi_new[lev][mfi].array();    
        const Box& box = mfi.tilebox();
        const int* lo      = box.loVect();
        const int* hi      = box.hiVect();
        const auto problo = Geom(lev).ProbLoArray();
        const auto dx     = Geom(lev).CellSizeArray();

        for(int n=0;n<Sborder.nComp();n++)
            for (int k=lo[2];k<=hi[2];k++)
                for (int j=lo[1];j<=hi[1];j++) 
                    for (int i=lo[0];i<=hi[0];i++)
                        {
                            Real x =  problo[0] + (i+0.5)*dx[0];
                            Real y =  problo[1] + (j+0.5)*dx[1];
                            Real z =  problo[2] + (k+0.5)*dx[2];

                            Real r2=x*x + y*y + z*z;

                            //std::cout << fab(i,j,k,n) << std::endl;
                            phi_new_fab(i,j,k,n)=
                            //fab(i,j,k,n);
        
        
                            //   std::pow( ( fab(i+1,j,k,n) - fab(i-1,j,k,n))/(dx[0]*2) ,2) +
                            // std::pow( ( fab(i,j+1,k,n) - fab(i,j-1,k,n))/(dx[1] *2),2) +
                            // std::pow( ( fab(i,j,k+1,n) - fab(i,j,k-1,n))/(dx[2]*2) ,2) ; 
                             

                             
                            //   ( fab(i+1,j,k,n) + fab(i-1,j,k,n) -2* fab(i,j,k,n) )/(dx[0]*dx[0] )+
                            // ( fab(i,j+1,k,n) + fab(i,j-1,k,n) -2* fab(i,j,k,n) )/(dx[1]*dx[1]) +
                            // ( fab(i,j,k+1,n) + fab(i,j,k-1,n) -2* fab(i,j,k,n) )/(dx[2]*dx[2]);  

                            ( -1./12. * fab(i-2,j,k,n) + 4./3 *fab(i-1,j,k,n) -5./2* fab(i,j,k,n) + 4/3.*fab(i+1,j,k,n)  -1./12*fab(i+2,j,k,n)  )/(dx[0]*dx[0] ) + 
                            ( -1./12. * fab(i,j-2,k,n) + 4./3 *fab(i,j-1,k,n) -5./2* fab(i,j,k,n) + 4/3.*fab(i,j+1,k,n)  -1./12*fab(i,j+2,k,n)  )/(dx[1]*dx[1] ) +
                            ( -1./12. * fab(i,j,k-2,n) + 4./3 *fab(i,j,k-1,n) -5./2* fab(i,j,k,n) + 4/3.*fab(i,j,k+1,n)  -1./12*fab(i,j,k+2,n)  )/(dx[2]*dx[2] );


                            





                        }
    }

    phi_new[lev].mult(dt_lev,num_grow);

    phi_new[lev].plus(phi_old[lev],0,phi_new[lev].nComp(),num_grow);

}


Real normCartesian( const MultiFab & phi_real ,  const Geometry & geom, int component)
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
        const int j=0;
        const int k=0;

#if AMREX_SPACEDIM == 1
            for (int i=lo[0];i<=hi[0];i++)
            {
             
                norm2+=(phi_real_box(i,j,0,component)*phi_real_box(i,j,0,component)
          ;
            }
    
#endif


#if AMREX_SPACEDIM == 2
        for (int j=lo[1];j<=hi[1];j++)
            for (int i=lo[0];i<=hi[0];i++)
            {
             
                norm2+=(phi_real_box(i,j,0,component)*phi_real_box(i,j,0,component)
          ;
            }
    
#endif




#if AMREX_SPACEDIM == 3
        for (int k=lo[2];k<=hi[2];k++)
            for (int j=lo[1];j<=hi[1];j++)
                for (int i=lo[0];i<=hi[0];i++)
                    {

                        norm2+=(phi_real_box(i,j,k,component)*phi_real_box(i,j,k,component) );
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


std::vector<Real> AmrCoreDiff::norm2( ) const
{
    auto & mfs = phi_new;


    int nComp=mfs[0].nComp();
    std::vector<Real> norms2(nComp,0);
    std::vector<Real> dV(max_level+1,1);
    for(int level=0;level<=max_level;level++)
    {
        const auto dx = Geom(level).CellSizeArray();
        dV[level]=1.;
        for (int d=0;d<amrex::SpaceDim;d++)
        {
            dV[level]*=dx[d];
        }

        
    }

    for (int level=0;level<=max_level;level++)
    {       
        for (int c=0;c<nComp;c++)
        {
            
            norms2[c]+=normCartesian(mfs[level],Geom(level),0) ;
        }
    }


    for (int level=max_level; level > 0    ;level--)
    {       
        BoxArray ba=mfs[level].boxArray();
        ba.coarsen(ref_ratio[level-1]);

        MultiFab overlappedRegion(ba, dmap[level], phi_new[level].nComp(), 0);

        
        overlappedRegion.copy(mfs[level-1]); // multifab containing the overlapped region
        for (int c=0;c<nComp;c++)
        {
            norms2[c]-=normCartesian(overlappedRegion,Geom(level-1),0);
        }

    }

    return norms2;

}