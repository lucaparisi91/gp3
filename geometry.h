#ifndef GEOMETRY_H
#define GEOMETRY_H


#include <algorithm>
#include<array>
#include "traits.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;
using Real = double;

class geometry
{
public:
	
	geometry( std::array<size_t,DIMENSIONS> shape_ ) : shape(shape_){  };

	std::array<size_t,DIMENSIONS> shape;
	std::array<Real,DIMENSIONS> lower_edges;
	std::array<Real,DIMENSIONS> higher_edges;
	int max_grid_size;
	
};

#endif