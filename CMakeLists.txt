cmake_minimum_required(VERSION 3.10)
# set the project name
project(gp)

include(cmake/init.cmake)

add_subdirectory(external/pybind11 EXCLUDE_FROM_ALL)
add_subdirectory(external/json EXCLUDE_FROM_ALL)
add_subdirectory(external/googletest EXCLUDE_FROM_ALL)
#add_subdirectory(external/pybind11_json EXCLUDE_FROM_ALL)

add_subdirectory(src)
add_subdirectory(tst)
add_subdirectory(pythonBindings)