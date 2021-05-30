
using Real = double;
#include <cmath>
#include <iostream>
#include "timers.h"
#include "mpiTools.h"








int main(int argc,char** argv)
{

    MPI_Init(&argc,&argv);

    
    int N=400;
    int nTrials = 10 ;
    if (argc >= 2)
    {
        N =std::atoi( argv[1]);
    }
    

    if (argc >= 3)
    {
        nTrials =std::atoi( argv[2]);
    }    

    int nGhosts=1;
    int nExtended = N + 2*nGhosts;

    




    Real * phi = new double[ nExtended*nExtended*nExtended];
    Real * phi2 = new double[ nExtended*nExtended*nExtended];


    Real L=10;


    int iLower=nGhosts;
    int iUpper=N-nGhosts;

    Real lowerEnd=-L/2.;
    Real upperEnd=L/2.;


    Real sum=0;
    Real deltax = (upperEnd - lowerEnd)/N;

    Real deltay=deltax;
    Real deltaz=deltax;

    Real alpha=1;

     for(int k=0;k<nExtended;k++)  
        for(int j=0;j<nExtended;j++)
            for(int i=0;i<nExtended;i++)
            {

                Real x=lowerEnd + (i - iLower + 0.5)*deltax;
                Real y=lowerEnd + (j - iLower + 0.5)*deltay;
                Real z=lowerEnd + (k - iLower + 0.5)*deltaz;
                auto r2=x*x + y*y + z*z;
                phi[k*nExtended*nExtended + j*nExtended + i ]=std::exp(-alpha*r2);   
                phi2[k*nExtended*nExtended + j*nExtended + i ]=0;

            }
    
    Real deltax2Inverse=1./(deltax*deltax);
    Real deltay2Inverse=1./(deltay*deltay);
    Real deltaz2Inverse=1./(deltaz*deltaz);

    std::cout << "nCellsPerGrid " << N <<  std::endl;



    for (int t=0;t<nTrials;t++)
    {
        START_TIMER("evaluate");

#pragma omp parallel for
    for(int k=iLower;k<iUpper;k++)  
        for(int j=iLower;j<iUpper;j++)
            for(int i=iLower;i<iUpper;i++)
            {
                phi2[k*nExtended*nExtended + j*nExtended + i ]=
                
                (- 2*phi[k*nExtended*nExtended + j*nExtended + i] 
                + phi[k*nExtended*nExtended + j*nExtended + i+1] 
                + phi[k*nExtended*nExtended + j*nExtended + i-1] )*deltax2Inverse + 
                (- 2*phi[k*nExtended*nExtended + j*nExtended + i] 
                + phi[k*nExtended*nExtended + (j+1)*nExtended + i] 
                + phi[k*nExtended*nExtended + (j-1)*nExtended + i] )*deltay2Inverse +
                (- 2*phi[k*nExtended*nExtended + j*nExtended + i] 
                + phi[(k+1)*nExtended*nExtended + j*nExtended + i] 
                + phi[(k-1)*nExtended*nExtended + j*nExtended + i] )*deltaz2Inverse

                ;

            }
        STOP_TIMER("evaluate");
    }

Real sums=0;

for(int k=0;k<nExtended*nExtended*nExtended;k++)  
{
    sums+=phi2[k]*phi2[k];
}

std::cout << sums << std::endl;

std::cout << timers::getInstance().report() << std::endl;

MPI_Finalize();
}