#ifndef MPI_TOOLS_H
#define MPI_TOOLS_H


#include "mpi.h"
#include <vector>
#include <memory>
#include <array>


using Real = double;


struct blockDecomposition
{

    public:
    using index_t = int;

    blockDecomposition( std::vector<index_t> nx ,std::vector<index_t> ny,std::vector<index_t> nz, std::array<index_t,2> nGhostsX,std::array<index_t,2> nGhostsY,std::array<index_t,2> nGhostsZ) : 
    _nx(nx) , _ny(ny) , _nz(nz),
    _nGhostsX(nGhostsX) , _nGhostsY(nGhostsY) , _nGhostsZ(nGhostsZ)
    {

    }  

    const auto & indices(int dim)
    {
        //assert(dim < 3);
        if (dim == 0 )
        {
            return _nx;
        }
        else if (dim == 1)
        {
            return _ny;
        }
        else
        {
            return _nz;
        }
    }

    int rank() ;


    const auto &  nx() const  {return _nx;}
    const auto &  ny() const  {return _ny;}
    const auto &  nz() const  {return _nz;}



    int nBlocks() const {
        return nBlocks(0)*nBlocks(1)*nBlocks(3);
    }

    int nBlocks( int i) const {

        switch(i) 
        {
            case 0 :
            {
                return _nx.size() - 1;
                break;
            }

            case 1 :
            {
                return _ny.size() - 1;
                break;
            }

            case 2 :
            {
                return _nz.size() - 1;
                break;
            }

        }
        return -1;
    }



    const auto & nGhostsX(){return _nGhostsX;}
    const auto & nGhostsY(){return _nGhostsY;}
    const auto & nGhostsZ(){return _nGhostsZ;}



    private:

    std::vector<index_t> _nx;
    std::vector<index_t> _ny;
    std::vector<index_t> _nz;


    std::array<index_t,2> _nGhostsX;
    std::array<index_t,2> _nGhostsY;
    std::array<index_t,2> _nGhostsZ;

    std::array<size_t,3> offset();

};



MPI_Datatype generateMPIDataTypeFace(int nx , int ny, int nz, int with, int dim);


blockDecomposition  generateBlockDecomposition(std::array<int,3> N , std::array<int,3> M ,
std::array<int,3> nGhosts={1,1,1});


class parallelComunication
{
    public:

    parallelComunication(MPI_Comm comm, std::shared_ptr<blockDecomposition> decomp );


    ~parallelComunication();


    private:
   

};

#endif