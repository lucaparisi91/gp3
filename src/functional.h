#ifndef GP_FUNCTIONAL_H
#define GP_FUNCTIONAL_H

#include "operators.h"

namespace gp
{

class functional
{

public:

    virtual void evaluate(
	wavefunction & waveNew, const wavefunction & waveOld, Real time=0)=0 ;
};

template<class laplacianOperator_t>
class trappedGPFunctional : public functional
{
    public:


    trappedGPFunctional(Real omega = 1, int nComp_=1) : prefactor(0.5 * omega*omega) , nComp(nComp_){};

	trappedGPFunctional(const json_t & j ) : prefactor(0.5*std::pow( j["omega"].get<Real>() , 2)    ) 
	  {

	}

    virtual void evaluate(
	wavefunction & waveNew, const  wavefunction & waveOld, Real time ) ;

	static std::string name() {return "harmonic";}

    private:

    Real prefactor;
	laplacianOperator_t _lap;
	int nSpecies;
	int offsetSpecies;
	int nComp;
};


template<class laplacianOperator_t>
class dropletFunctional : public functional
{
    public:

    dropletFunctional() {};
	dropletFunctional(const json_t & j ){}

    virtual void evaluate(
	wavefunction & waveNew, const  wavefunction & waveOld, Real time ) ;

	static std::string name() {return "droplet";}

    private:

	laplacianOperator_t _lap;

};

void setGaussian(
	wavefunction & waveNew, Real alpha, int c);

}


#endif