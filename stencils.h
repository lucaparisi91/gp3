
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
	static inline Real x(  T & phi, int i , int j, const Real  dx)
	{
		return  ((phi(i+1,j,0) - phi(i-1,j,0) )/(2*dx) ); 
	};

	public:
	template<class T>
	static inline Real y(  T & phi, int i , int j, const Real  dx)
	{
		return  (phi(i,j+1,0) - phi(i,j-1,0) )/(2*dx) ; 
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
	static inline Real x(  T & phi, int i , int j, const Real  dx)
	{
		return (-3./2*phi(i,j,0)	+2*phi(i+1,j,0) - 1./2 * phi(i+2,j,0))/dx ;
	};

};



template<>
class laplacian<2,CENTRAL>
{
public:
	template<class T>
	static inline Real x(  T & phi, int i , int j, const Real  dx)
	{
		return  (phi(i-1,j,0) + phi(i+1,j,0) - 2*phi(i,j,0) )/(dx*dx) ; 
	};
	template<class T>
	static inline Real y(  T & phi, int i , int j, const Real dx)
	{
		return  (phi(i,j-1,0) + phi(i,j+1,0) - 2*phi(i,j,0) )/(dx*dx) ; 
	};

	template<class T>
	static inline Real call(  T & phi, int i , int j, const Real *dx)
	{
		return  x(phi,i,j,dx[0]) + y(phi,i,j,dx[1]); 
	};

};

template<>
class laplacian<2,FORWARD>
{
public:
	template<class T>
	static inline Real x(  T & phi, int i , int j, const Real  dx)
	{
		return ( 2*phi(i,j,0) -5*phi(i+1,j,0) + 4*phi(i+2,j,0)	-1 * phi(i+3,j,0) )/(dx*dx);	 
	};

};

template<int order,int dir>
class laplacianCylindricalSymm
{
public:
	template<class T>
	static inline Real x(  T & phi, int i , int j, const Real dx, const Real x)
	{
		return 
		laplacian<order,dir>::x(phi,i,j,dx) + gradient<order,dir>::x(phi,i,j,dx)*1./x
		;
	};

	template<class T>
	static inline Real y(  T & phi, int i , int j, const Real dx)
	{
		return laplacian<order,dir>::y(phi,i,j,dx);
	};



};


