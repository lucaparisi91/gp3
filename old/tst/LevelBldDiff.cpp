
#include <AMReX_LevelBld.H>
#include "AmrLevelDiff.h"

using namespace amrex;

class LevelBldDiff
    :
    public LevelBld
{
    virtual void variableSetUp () override;
    virtual void variableCleanUp () override;
    virtual AmrLevel *operator() () override;
    virtual AmrLevel *operator() (Amr&            papa,
                                  int             lev,
                                  const Geometry& level_geom,
                                  const BoxArray& ba,
				  const DistributionMapping& dm,
                                  Real            time) override;
};

LevelBldDiff Diff_bld;

LevelBld*
getLevelBld ()
{
    return &Diff_bld;
}

void
LevelBldDiff::variableSetUp ()
{
    AmrLevelDiff::variableSetUp();
}

void
LevelBldDiff::variableCleanUp ()
{
    AmrLevelDiff::variableCleanUp();
}

AmrLevel*
LevelBldDiff::operator() ()
{
    return new AmrLevelDiff;
}

AmrLevel*
LevelBldDiff::operator() (Amr&            papa,
	   	         int             lev,
                         const Geometry& level_geom,
                         const BoxArray& ba,
                         const DistributionMapping& dm,
                         Real            time)
{
    return new AmrLevelDiff(papa, lev, level_geom, ba, dm, time);
}
