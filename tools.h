#ifndef TOOLS_H
#define TOOLS_H
#include <AMReX_PlotFileUtil.H>
using namespace amrex;

Real norm( const MultiFab & phi_real , const MultiFab & phi_imag,  const Geometry & geom, int component=0);



class initializer
{
    bool isInitialized;
    static initializer *s_instance;
    initializer() : isInitialized(false ) {}
  

  public:
    void init()
    {
        if (!isInitialized)
        {
            amrex::Initialize(MPI_COMM_WORLD);
            isInitialized = true;
        }
    }



    static initializer  & instance()
    {
        if (!s_instance)
          s_instance = new initializer;
        return *s_instance;
    }
};

// Allocating and initializing GlobalClass's
// static data member.  The pointer is being
// allocated - not the object inself.

#endif