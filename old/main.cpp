#include "run.h"
using Real = double;
#include <AMReX_PlotFileUtil.H>

int main(int argc,char** argv)
{
	
  run();
  amrex::Finalize();

}
