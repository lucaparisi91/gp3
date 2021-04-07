#include "timers.h"
#include <iostream>
#include <sstream> 
#include <mpi.h>
timers & timers::getInstance()
{
    if (_singleton == nullptr )
    {
        _singleton= new timers();
    }

    return (*_singleton);
}

timers*  timers::_singleton = nullptr;




std::string timers::report() const
{
    std::stringstream ss;
    int rank;
    int nProcs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD,&nProcs);
     ss << "Rank " << rank << std::endl;
     ss << "------------------------" << std::endl;
    for (auto it : _timers)

    {
       
        ss << "Timer " <<it.first <<": " << (it.second)->timeElapsed() << std::endl; 
    }
    
    std::string local_report=ss.str();
    
   

    int nCharactersLocal=local_report.size();
    
    int nCharacters[nProcs];
    MPI_Allgather(
        &nCharactersLocal,
        1,
        MPI_INT,
        & nCharacters[0],
        1,
        MPI_INT,
        MPI_COMM_WORLD
    );

    int displs[nProcs];
    displs[0]=0;

    for(int i=1;i<nProcs;i++)
    {
        displs[i]=displs[i-1] + nCharacters[i-1];
    }

    int nCharactersGlobal=0;

    for(int i=0;i<nProcs;i++)
    {
        nCharactersGlobal+=nCharacters[i];
    }


    char globalString[nCharactersGlobal+1];

    MPI_Allgatherv(
    & (local_report[0]),
    nCharactersLocal,
    MPI_CHAR,
    globalString,
    &(nCharacters[0]),
    displs,
    MPI_CHAR,
    MPI_COMM_WORLD);

    globalString[nCharactersGlobal]='\0';

    std::string globalReport(globalString);


    return globalReport;
}