#ifndef INITIALIZER_H
#define INITIALIZER_H

#include"traits.h"
#include "functionalFactory.h"
#include <AMReX_PlotFileUtil.H>
#include "operatorsFactory.h"


class initializer
{
    bool isInitialized;
    static initializer *s_instance;
    initializer();
    
    functionalFactory func;
    opFactory operatorFac;
  public:

    auto & getFunctionalFactory() {return func;}

    auto & getOperatorsFactory() {return operatorFac;}

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
// allocated - not the object inself.s


#endif