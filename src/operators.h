
#ifndef OPERATORS_H
#define OPERATORS_H
#include "wavefunction.h"

namespace gp
{

namespace operators
{


struct operatorBase {};

class reductionOperator : public operatorBase{};


template<int order,int dim>
struct laplacian : public operatorBase
{
    laplacian() {
        static_assert(order == 1,"Only 1st order laplacian is supported");
        static_assert(dim == 3, "Only 3D laplacians are supported");
    }

    void setGeometry( const geometry & geom)
    {
       
        for (int d=0;d<DIMENSIONS;d++)
        {
            laplacianMultiplicationFactor[d]=1./(geom.CellSize()[d]*geom.CellSize()[d]);
        }

    }


    Real operator()(int i, int j , int k, int c, const gp::const_array_t & phi) const
    {   
        if constexpr (dim == 3 and order==1 )
        {
            return ( phi(i+1,j,k,c) -2*phi(i,j,k,c) + phi(i-1,j,k,c) )*laplacianMultiplicationFactor[0]+
                    ( phi(i,j+1,k,c) -2*phi(i,j,k,c) + phi(i,j-1,k,c) )*laplacianMultiplicationFactor[1]+ 
                    ( phi(i,j,k+1,c) -2*phi(i,j,k,c) + phi(i,j,k-1,c) )*laplacianMultiplicationFactor[2];
        }

    }


    Real operator()(int i,int j, int k,int c, wavefunctionRegion & wave) const
    {
        const auto & phi = wave.getPhi<gp::array_t>();
        const auto & geom = wave.getGeometry();

        return (*this)(i,j,k,c,phi);
        
    }
    private:

     std::array<Real,3 > laplacianMultiplicationFactor;

};

struct setPhi : public operatorBase
{

    template<class op_t,class ...Args>
    void operator()(int i ,int j, int k,int c, wavefunctionRegion &  wave,op_t && currentOp,Args && ... args )
    {
        const auto & phi = wave.getPhi<gp::array_t>();


        phi(i,j,k,c)=currentOp(i,j,k,c,wave,std::forward<Args>(args)...);


    };


    template<class op_t,class ...Args>
    void operator()(int i ,int j, int k,int c,wavefunctionRegion &  waveNew, wavefunctionRegion &  waveOld,op_t && currentOp,Args && ... args )
    {
        const auto & phi = waveNew.getPhi<gp::array_t>();

        phi(i,j,k,c)=currentOp(i,j,k,c,waveOld,std::forward<Args>(args)...); 

    };

};



template<class T1,class T2>
class sumOperators : public operatorBase
{
    public:
    sumOperators(T1 left,T2 right) : _left(left),_right(right) {}


    template<class ... Args>
    auto operator()(Args&& ... args) {

        
        return _left(std::forward<Args>(args)...)   + _right(std::forward<Args>(args)...)  ; 
        }
    
    private:

    T1 _left;
    T2 _right;
}; 




template<class T1,class T2>
class productOperators : public operatorBase
{
    public:
    productOperators(T1 left,T2 right) : _left(left),_right(right) {}

    
    template<class ... Args>
    auto operator()(Args&& ... args) {

        
        return _left(std::forward<Args>(args)...)   * _right(std::forward<Args>(args)...)  ; 
        }
    
    private:

    T1 _left;
    T2 _right;
}; 

template<class T1,class T2>
class productScalarOperator : public operatorBase
{
    public:
    productScalarOperator(T1 left,T2 right) : _left(left),_right(right) {}

    
    template<class ... Args>
    auto operator()(Args&& ... args) {
        return _left   * _right(std::forward<Args>(args)...)  ; 
        }
    
    private:

    T1 _left;
    T2 _right;
}; 



template<class T1,class T2, std::enable_if_t< std::is_base_of<operatorBase,T1>::value  >* =nullptr   ,std::enable_if_t< std::is_base_of<operatorBase,T2>::value  >* =nullptr  >
sumOperators<T1,T2> operator+(T1 left, T2 right)
{
    return sumOperators<T1,T2>(left,right);
}

template<class T1,class T2, std::enable_if_t< !std::is_base_of<operatorBase,T1>::value  >* =nullptr   ,std::enable_if_t< std::is_base_of<operatorBase,T2>::value  >* =nullptr  >
auto operator*(T1 left, T2 right)
{
    return productScalarOperator<T1,T2>(left,right);
}


};


};

#endif


