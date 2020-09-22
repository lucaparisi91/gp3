#include "geometry.h"
#include "run.h"
#include "evaluate.h"
#include "model.h"
#include "tools.h"

py::array_t<double, py::array::f_style | py::array::forcecast> evaluatePython(py::array_t<double>  real, py::array_t<double>  imag ,  const geometry & geom )
{
		 initializer::instance().init();

         int order = 4;

		 model m(geom, order , 1);
         model m2(geom, order , 1);
         
		 m.fill(real,imag);
         
         evaluate(m2.real(), m2.imag() ,
            m.real() , m.imag() , 
            0 , m.getGeometry() , m.realLaplacian() , m.imagLaplacian() 
            );
     

         py::array_t<double> values( geom.size() );
         size_t t =0;
         double * p = (double * )(values.request().ptr );

         LOOP(m2.real() , m2.getGeometry() )

         *(p + t ) = data(i,j,k,0);
         t++;
         ENDLOOP
         
        values.resize( geom.shape );

         
        return values;

		 
}


PYBIND11_MODULE(gp_c,m) {
m.doc() = "GP simulation python binding"; // optional module docstring
py::class_<geometry>(m, "geometry")
    .def(py::init<std::array<size_t,DIMENSIONS> >())
    .def_readwrite("shape", &geometry::shape)
    .def_readwrite("lower_edges", &geometry::lower_edges)
    .def_readwrite("higher_edges", &geometry::higher_edges)
    .def_readwrite("max_grid_size", &geometry::max_grid_size)
    .def_readwrite("symmetry", &geometry::symmetry);
m.def("run", &run, "Run the simulation");
m.def("evaluate",&evaluatePython, "Evaluate the rhs." );

}


