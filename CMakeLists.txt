cmake_minimum_required(VERSION 3.10)
# set the project name
project(gp)

include(cmake/init.cmake)

add_subdirectory(external/pybind11)
add_subdirectory(external/amrex)
add_subdirectory(external/json)
add_subdirectory(external/pybind11_json)
add_subdirectory(external/googletest)

add_subdirectory(src)
add_subdirectory(tst)
add_subdirectory(pythonBindings)