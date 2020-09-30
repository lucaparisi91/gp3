
enum {CENTRAL = 0 , FORWARD = 1, BACKWARD= 2} directionDeriv;

template<int i , int dir=CENTRAL>
class laplacian;

template<int i , int dir=CENTRAL>
class gradient;

template<>
class gradient<2,CENTRAL>
{
	public:
	template<class T>
	static inline Real x(  T & phi, int i , int j, int k,const Real  dx)
	{
		return  ((phi(i+1,j,k) - phi(i-1,j,k) )/(2*dx) ); 
	};

	public:
	template<class T>
	static inline Real y(  T & phi, int i , int j,int k, const Real  dx)
	{
		return  (phi(i,j+1,k) - phi(i,j-1,k) )/(2*dx) ; 
	};

};

template<>
class gradient<4,CENTRAL>
{
	public:
	template<class T>
	static inline Real x(  T & phi, int i , int j, const Real  dx)
	{
		return  ( 
			1./12 * phi(i-2,j,0)  -2./3* phi(i-1,j,0) 
			-1./12 * phi(i+2,j,0)  +2./3* phi(i+1,j,0) 
			)/(dx) ; 
	};
};


template<>
class gradient<2,FORWARD>
{
	public:
	template<class T>
	static inline Real x(  T & phi, int i , int j, int k,const Real  dx)
	{
		return (-3./2*phi(i,j,k)	+2*phi(i+1,j,k) - 1./2 * phi(i+2,j,k))/dx ;
	};



};



template<>
class laplacian<2,CENTRAL>
{
public:
	template<class T>
	static inline Real x(  T & phi, int i , int j, int k , const Real  dx)
	{
		return  (phi(i-1,j,k) + phi(i+1,j,k) - 2*phi(i,j,k) )/(dx*dx) ; 
	};
	template<class T>
	static inline Real y(  T & phi, int i , int j, int k, const Real dx)
	{
		return  (phi(i,j-1,k) + phi(i,j+1,k) - 2*phi(i,j,k) )/(dx*dx) ; 
	};

	template<class T>
	static inline Real z(  T & phi, int i , int j, int k, const Real dx)
	{
		return  (phi(i,j,k-1) + phi(i,j+1,k+1) - 2*phi(i,j,k) )/(dx*dx) ; 
	};


	template<class T>
	static inline Real call(  T & phi, int i , int j, int k, const Real *dx)
	{

#if AMREX_SPACEDIM == 2
		return  x(phi,i,j,k,dx[0]) + y(phi,i,j,k,dx[1]); 
#endif
#if AMREX_SPACEDIM == 1
		return  x(phi,i,j,k,dx[0]); 
#endif

#if AMREX_SPACEDIM == 3
		return  x(phi,i,j,k,dx[0]) + y(phi,i,j,k,dx[1]) +  y(phi,i,j,k,dx[2]); 
#endif

	};

};

template<>
class laplacian<2,FORWARD>
{
public:
	template<class T>
	static inline Real x(  T & phi, int i , int j, int k,const Real  dx)
	{
		return ( 2*phi(i,j,k) -5*phi(i+1,j,k) + 4*phi(i+2,j,k)	-1 * phi(i+3,j,k) )/(dx*dx);	 
	};

};

template<int order,int dir>
class laplacianCylindricalSymm
{
public:
	template<class T>
	static inline Real x(  T & phi, int i , int j, int k, const Real dx, const Real x)
	{
		return 
		laplacian<order,dir>::x(phi,i,j,k,dx) + gradient<order,dir>::x(phi,i,j,k,dx)*1./x
		;
	};

	template<class T>
	static inline Real y(  T & phi, int i , int j, int k, const Real dx)
	{
		return laplacian<order,dir>::y(phi,i,j,k,dx);
	};



};


template<int order,int dir>
class laplacianSphericalSymm
{
public:

	template<class T>
	static inline Real call(  T & phi, int i , int j, int k, const Real *dx, const Real x)
	{
		#if AMREX_SPACEDIM == 1
		return laplacian<order,dir>::x(phi,i,j,k,dx[0]) + 2*gradient<order,dir>::x(phi,i,j,k,dx[0])*1./x
		;
		#endif

	};



};






