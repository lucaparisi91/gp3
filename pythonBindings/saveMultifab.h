#include <pybind11/pybind11.h>
#include "pybind11_json/pybind11_json.hpp"
#include <pybind11/numpy.h>
#include <nlohmann/json.hpp>
#include <vector>

namespace py = pybind11;
using json_t = nlohmann::json ;
using Real = double;


void saveMultifab(std::vector< py::array_t<std::complex<Real> > > initialConditions  , const json_t & settings   );