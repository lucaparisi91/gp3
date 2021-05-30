// one dimensional decomposition in uniform slabs

#include "mpiTools.h"

blockDecomposition generateBlockDecomposition(std::array<int,3> N , std::array<int,3> M ,
std::array<int,3> nGhosts)
{
    constexpr int dim=3;
    std::array<std::vector<int>,dim > partitions;

   
    

    for (int d=0;d<dim;d++)
    {
        auto & nx = partitions[d];

        int deltan=M[d]/N[d];
        int remainder=M[d]%N[d];
        std::vector<size_t> spacing(N[d] ,deltan);

        for(int n=N[d]  -remainder;n<N[d] ;n++)
        {
            spacing[n]+=1;
        };

        nx.resize(N[d]+1,0);

        for(int i=1;i<N[d] + 1;i++)
        {
            nx[i]=nx[i-1] + spacing[i-1];
        }

    }

    return blockDecomposition(partitions[0],partitions[1],partitions[2] , {nGhosts[0],nGhosts[0] },{nGhosts[1],nGhosts[1] },{nGhosts[2],nGhosts[2] } );

    
}


MPI_Datatype generateMPIDataTypeFace(int nx , int ny, int nz, int width, int dim)
{
    MPI_Datatype faceType;

    if (dim == 0 )

    {
        int err=MPI_Type_vector(ny*nz,width,nx,MPI_DOUBLE, & faceType );
    }
    else if (dim == 1)
    {
        int err=MPI_Type_vector(nz,nx*width,ny,MPI_DOUBLE, & faceType );
    }
    else if (dim == 2)
    {
        int err=MPI_Type_vector(width,nx*ny,0,MPI_DOUBLE, & faceType );
    }
    MPI_Type_commit(&faceType);

    return faceType;
}

MPI_Datatype delteMPIDataTypeFace(MPI_Datatype & faceDataType)
{
    MPI_Type_free(&faceDataType);
    return faceDataType;
}