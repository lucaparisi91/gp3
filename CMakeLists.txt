cmake_minimum_required(VERSION 3.10)
# set the project name
project(gp)
set(DIM 3)


set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_INSTALL_PREFIX /home/luca/source/GP3/install)

list(APPEND CMAKE_PREFIX_PATH CMAKE_PREFIX_PATH /home/luca/source/GP3/external/arrow_cpp)

add_subdirectory(external/pybind11)
add_subdirectory(external/amrex)
add_subdirectory(external/json)
add_subdirectory(external/pybind11_json)
add_subdirectory(external/googletest)

find_package(Arrow REQUIRED)

link_directories(  /home/luca/source/GP3/external/arrow_cpp/lib)
include_directories( /home/luca/source/GP3/external/arrow_cpp/include)


#include_directories(AMReX/include)
#link_directories(AMReX/lib)
find_package(MPI REQUIRED)

find_package(PythonLibs)
include_directories(${PYTHON_INCLUDE_DIRS})


set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}   )
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS}  )
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wfatal-errors -std=gnu++17")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wfatal-errors -std=gnu++17")

set(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY TRUE)
message(${CMAKE_SOURCE_DIR})
include_directories(${MPI_INCLUDE_PATH} external/json/single_include external/pybind11_json/include)




#set( CMAKE_PREFIX_PATH AMReX)
set(GP_FILES run.cpp stepper.cpp geometry.cpp evaluate.cpp model.cpp tools.cpp operators.cpp
functional.cpp initializer.cpp plotFile.cpp gpDropletFunctional3D.cpp)


add_library(gp_lib ${GP_FILES} )
target_link_libraries(gp_lib PRIVATE amrex ${MPI_LIBRARIES} gfortran arrow parquet  )

target_include_directories(gp_lib PUBLIC
                           external/pybind11/include
                           external/json/include
                           external/pybind11_json/include
			   )

set(pythonBindingTarget gp${DIM}D_c)

pybind11_add_module(${pythonBindingTarget} gp.cpp )
#find_package(AMReX)
target_link_libraries( ${pythonBindingTarget} PRIVATE amrex ${MPI_LIBRARIES} gfortran gp_lib arrow )
#add_executable(main main.cpp run.cpp)
#target_link_libraries( main PUBLIC amrex ${MPI_LIBRARIES} gfortran pthread)

set (TEST_SOURCES tst/model-test.cpp tst/test-multilevel.cpp)
add_executable(tst tst/main.cpp ${TEST_SOURCES})
target_link_libraries(tst PRIVATE  gtest  ${MPI_LIBRARIES} amrex arrow parquet)
target_link_libraries(tst PUBLIC  gp_lib)

target_include_directories(tst PUBLIC
                           "${PROJECT_BINARY_DIR}/.."
                           )