#include "geometry.h"
#include "run.h"
#include "evaluate.h"
#include "model.h"
#include "tools.h"
#include "pybind11_json/pybind11_json.hpp"
#include "initializer.h"
#include "functional.h"


auto  toPyArray(MultiFab & real, MultiFab & imag, Geometry & geom)
{
    auto shape = geom.Domain().size();
    size_t size=1;

    for (int d=0;d<AMREX_SPACEDIM;d++)
    {
        size*=shape[d];
    }

    std::cout<< size << std::endl;

    using return_type =  py::array_t<std::complex<double> , py::array::f_style | py::array::forcecast  > ;


    return_type values( size );

    values.resize(  {AMREX_D_DECL( shape[0] , shape[1] , shape[2]  ) });

    std::complex<double> * p = (std::complex<double> * )(values.request().ptr );
    
    size_t t=0;
    LOOP( real , geom )

         *(p + t ) = data (  i, j, k,0) + 1i*0.0;
         t++;
   ENDLOOP

   t=0;

    LOOP( imag , geom )
         *(p + t ) += data(i,j,k,0) *1i;
         t++;
   ENDLOOP

   return values;


}



auto   evaluatePython( py::array_t<std::complex<Real> > initialCondition , json_t & settings  )
{
	initializer::instance().init();

    auto [ box, geom , dm, low_bc, high_bc] = createGeometry(settings["geometry"]);
    
    int Ncomp = 1;
    int order = settings["functional"]["laplacian"]["order"];
    int Nghost = 1;

    MultiFab phi_real_old(box, dm, Ncomp, Nghost);
    MultiFab phi_imag_old(box, dm, Ncomp, Nghost);

    MultiFab phi_real_new(box, dm, Ncomp, Nghost);
    MultiFab phi_imag_new(box, dm, Ncomp, Nghost);



    phi_real_old=0.;
    phi_imag_old=0.;

    phi_real_new=0.;
    phi_imag_new=0.;

    fill(phi_real_old, phi_imag_old, initialCondition , geom);

    auto func = initializer::instance().getFunctionalFactory().create(settings["functional"]);
    
    
    func->define(geom,box,dm,low_bc,high_bc);
    phi_real_old.setDomainBndry(0, 0, 1, geom)   ;
    phi_imag_old.setDomainBndry(0, 0, 1, geom)   ;


    auto bc = toMultiFabBC(low_bc,high_bc);


    func->evaluate(phi_real_new,phi_imag_new,phi_real_old,phi_imag_old,0);

    delete func;

    return toPyArray(phi_real_new, phi_imag_new,geom);
    
} 

#define _MODULE_NAME_(n) gp##n ## D_c
#define MODULE_NAME(n) _MODULE_NAME_(n)

PYBIND11_MODULE( MODULE_NAME(AMREX_SPACEDIM)   ,m) {
m.doc() = "GP simulation python binding"; // optional module docstring
/* py::class_<geometry>(m, "geometry")
    .def(py::init<std::array<size_t,AMREX_SPACEDIM> >())
    .def_readwrite("shape", &geometry::shape)
    .def_readwrite("lower_edges", &geometry::lower_edges)
    .def_readwrite("higher_edges", &geometry::higher_edges)
    .def_readwrite("max_grid_size", &geometry::max_grid_size)
    .def_readwrite("symmetry", &geometry::symmetry);
 m.def("run", &run, "Run the simulation");
 */
m.def("evaluate",&evaluatePython, "Evaluate the rhs." );
m.def("run", & run , "Run the GP simulation");
}


