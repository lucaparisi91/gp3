#include "gtest/gtest.h"
#include "mpiTools.h"


TEST(mpi, decomposition)
{
    {
    const int M = 100;
    const int N = 2;

    auto decomp = generateBlockDecomposition({N,1,1},{M,M,M});

    ASSERT_EQ(decomp.nx().size() , 3);
    ASSERT_EQ(decomp.nx()[0] , 0 );
    ASSERT_EQ(decomp.nx()[1] , 50 );
    ASSERT_EQ(decomp.nx()[2] , 100 );
    ASSERT_EQ(decomp.ny()[0] , 0 );
    ASSERT_EQ(decomp.ny()[1] , 100 );



    }

    {
    const int M = 100;
    const int N = 3;

    auto decomp = generateBlockDecomposition({N,1,1},{M,M,M});

    ASSERT_EQ(decomp.nx().size() , N + 1);
    ASSERT_EQ(decomp.nx()[0] , 0 );
    ASSERT_EQ(decomp.nx()[1] , 33 );
    ASSERT_EQ(decomp.nx()[2] , 66 );
    ASSERT_EQ(decomp.nx()[3] , 100 );
    ASSERT_EQ(decomp.ny()[0] , 0 );
    ASSERT_EQ(decomp.ny()[1] , 100 );
    ASSERT_EQ(decomp.nz()[0] , 0 );
    ASSERT_EQ(decomp.nz()[1] , 100 );
    
    }

    {

    auto decomp = generateBlockDecomposition({2,2,1},{512,512,512});


    ASSERT_EQ(decomp.nx()[0] , 0 );
    ASSERT_EQ(decomp.nx()[1] ,  256);
    ASSERT_EQ(decomp.nx()[2] , 512 );
    ASSERT_EQ(decomp.ny()[0] , 0 );
    ASSERT_EQ(decomp.ny()[1] , 256 );
    ASSERT_EQ(decomp.ny()[2] , 512 );
    ASSERT_EQ(decomp.nz()[0] , 0 );
    ASSERT_EQ(decomp.nz()[1] , 512 );
    
    }


}

TEST(mpi, initialization)
{
   /*  {

    std::array<Real,3> lBox{1,1,1};
    std::array<Real,3> leftEdge{0.5,-0.5,-0.5};

    const int nGhosts=1;
    std::array<Real,3> lBox;
    std::array<int,3> deltax{lBox[0]/M,lBox[1]/M,lBox[2]/M};    

    auto decomp = generateSlabDecomposition( M,nGhosts,MPI_COMM_WORLD);


    int nExtended = M + 2*nGhosts;

    const auto * data = new double[nExtended*nExtended*nExtended];

    std::array<int,3> lower=decomp.lowerIndexValid();
    std::array<int,3> upper=decomp.upperIndexValid();
    std::array<int,3> offset=decomp.offset();

    for (int i=lower[0];i<=upper[0];i++)
        for (int j=lower[1];j<=upper[1];j++)
            for (int k=lower[2];k<=upper[2];k++)
             {
                size_t index = k*M[0]*M[1] + j*M[0] + i;
                
                Real x = lowerEdge[0] + (i + offset[0]  - lower[0])*deltax[0];
                Real y = lowerEdge[1] + (j + offset[1] - lower[1])*deltax[1];
                Real z = lowerEdge[2] + (k + offset[2] - lower[2])*deltax[2];
             }

    decomp.sendReceiveFace(0,1,0,1);

 */


}

