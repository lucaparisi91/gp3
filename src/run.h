#ifndef RUN_H
#define RUN_H
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
namespace py = pybind11;
#include "geometry.h"
#include <nlohmann/json.hpp>

using json_t = nlohmann::json ;

void run(py::array_t<std::complex<Real> > initialConditions , const json_t & settings   );



#endif