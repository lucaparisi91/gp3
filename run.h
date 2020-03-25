#ifndef RUN_H
#define RUN_H
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
namespace py = pybind11;
#include "geometry.h"


void run(py::array_t<double> initialCondition_real,py::array_t<double> initialCondition_imag,const geometry & geomInfo);


#endif
