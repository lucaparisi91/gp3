

namespace gp
{

namespace operators
{
#include "wavefunction.h"


struct operatorBase {};

class reductionOperator : public operatorBase{};



template<int order,int dim>
struct laplacian : public operatorBase
{
    laplacian() {
        static_assert(order == 1,"Only 1st order laplacian is supported");
        static_assert(dim == 3, "Only 3D laplacians are supported");

    }
    Real operator()(int i,int j, int k,const wavefunctionRegion & wave,int c) const
    {
        const auto & phi = wave.getPhiOld<gp::array_t>();
        const auto & geom = wave.getGeometry();


        if constexpr (dim == 3 and order==1 )
        {
            return ( phi(i+1,j,k,c) -2*phi(i,j,k,c) + phi(i-1,j,k,c) )/(geom.CellSize()[0]*geom.CellSize()[0]) +
                    ( phi(i,j+1,k,c) -2*phi(i,j,k,c) + phi(i,j-1,k,c) )/(geom.CellSize()[1]*geom.CellSize()[1]) + 
                    ( phi(i,j,k+1,c) -2*phi(i,j,k,c) + phi(i,j,k-1,c) )/(geom.CellSize()[2]*geom.CellSize()[2]) ;
        }

    }
};

struct setPhiNew : public operatorBase
{

    template<class op_t,class ...Args>
    void operator()(int i ,int j, int k,wavefunctionRegion &  wave,op_t && currentOp,Args && ... args )
    {
        const auto & phi = wave.getPhiNew<gp::array_t>();

        phi(i,j,k,0)=currentOp(i,j,k,wave,0,std::forward<Args>(args)...); 


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