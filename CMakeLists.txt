cmake_minimum_required(VERSION 3.10)
# set the project name
project(gp)
add_subdirectory(external/pybind11)
include_directories(AMReX/include)
link_directories(AMReX/lib)
find_package(MPI REQUIRED)

find_package(PythonLibs)
include_directories(${PYTHON_INCLUDE_DIRS})

set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
message(${CMAKE_SOURCE_DIR})
include_directories(${MPI_INCLUDE_PATH})

set( CMAKE_PREFIX_PATH AMReX)
find_package(AMReX)
pybind11_add_module(gp_c gp.cpp run.cpp stepper.cpp evaluate.cpp)
target_link_libraries(gp_c PRIVATE amrex ${MPI_LIBRARIES} gfortran )
#add_executable(main main.cpp run.cpp)
#target_link_libraries( main PUBLIC amrex ${MPI_LIBRARIES} gfortran pthread)


install(TARGETS gp_c DESTINATION /home/luca/source/GP3/python)

