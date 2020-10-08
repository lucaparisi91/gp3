#include "functional.h"


class gpDroplet3D : public functional
{
    public:

    gpDroplet3D()  {};

	gpDroplet3D(const json_t & j ) :   functional::functional(j)
	{

	}


    virtual void evaluate(
	MultiFab & state_new_real, MultiFab & state_new_imag,
	MultiFab & state_old_real,  MultiFab & state_old_imag, Real time =0 ) ;

	static std::string name() {return "gpDroplet3D";}

    private:

};
