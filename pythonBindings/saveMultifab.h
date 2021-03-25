#include <pybind11/pybind11.h>
#include "pybind11_json/pybind11_json.hpp"
#include <pybind11/numpy.h>
#include <nlohmann/json.hpp>
#include <vector>
#include <pybind11/stl.h>
#include <pybind11/complex.h>


namespace py = pybind11;
using json_t = nlohmann::json ;
using Real = double;

void saveMultifab( py::list initialConditions   , const json_t & settings   );

std::vector< std::vector<std::complex<double> > > readMultifab( const json_t & settings   );

