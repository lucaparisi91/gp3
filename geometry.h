#ifndef GEOMETRY_H
#define GEOMETRY_H


#include <algorithm>
#include<array>
#include "src/traits.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;
using Real = double;
using RealArray = std::array<Real,AMREX_SPACEDIM>;
using IntArray = std::array<size_t,AMREX_SPACEDIM>;


class geometry
{
public:
	geometry( std::array<size_t,AMREX_SPACEDIM> shape_ ) : shape(shape_),symmetry("none"){  };
	size_t size() const;
	std::array<size_t,AMREX_SPACEDIM> shape;
	std::array<Real,AMREX_SPACEDIM> lower_edges;
	std::array<Real,AMREX_SPACEDIM> higher_edges;
	int max_grid_size;
	std::string symmetry;
	private:

};



#endif