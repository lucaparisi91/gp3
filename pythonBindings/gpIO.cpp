#include "saveMultifab.h"





PYBIND11_MODULE(gpIO_c, m) {
    m.doc() = "IO interface for GP data"; // optional module docstring

    m.def("saveMultifab", &saveMultifab, "arguments: json file containting geometry, boxes and a list of numpy arrays");
    m.def("readMultifab", &readMultifab, "arguments: json file containting geometry, boxes and name");

}