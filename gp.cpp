#include "geometry.h"
#include "run.h"


PYBIND11_MODULE(gp_c,m) {
m.doc() = "GP simulation python binding"; // optional module docstring
py::class_<geometry>(m, "geometry")
    .def(py::init<std::array<size_t,DIMENSIONS> >())
    .def_readwrite("shape", &geometry::shape)
    .def_readwrite("lower_edges", &geometry::lower_edges)
    .def_readwrite("higher_edges", &geometry::higher_edges)
    .def_readwrite("max_grid_size", &geometry::max_grid_size); 
m.def("run", &run, "Run the simulation");
}

