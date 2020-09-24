cmake_minimum_required(VERSION 3.10)
# set the project name
project(gp)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_INSTALL_PREFIX /home/luca/source/GP3/install)
set(DIM 3)
set(CMAKE_CXX_STANDARD 17)


add_subdirectory(external/pybind11)
add_subdirectory(external/amrex)
add_subdirectory(external/json)
add_subdirectory(external/pybind11_json)
add_subdirectory(external/googletest)


	           




#include_directories(AMReX/include)
#link_directories(AMReX/lib)
find_package(MPI REQUIRED)

find_package(PythonLibs)
include_directories(${PYTHON_INCLUDE_DIRS})


set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS} -Wfatal-errors)
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS}  )
message(${CMAKE_SOURCE_DIR})
include_directories(${MPI_INCLUDE_PATH} external/json/single_include external/pybind11_json/include)

#set( CMAKE_PREFIX_PATH AMReX)
set(GP_FILES run.cpp stepper.cpp geometry.cpp evaluate.cpp model.cpp tools.cpp operators.cpp
functional.cpp )

add_library(gp_lib ${GP_FILES} )
target_link_libraries(gp_lib PRIVATE amrex ${MPI_LIBRARIES} gfortran  )

target_include_directories(gp_lib PUBLIC
                           external/pybind11/include
                           external/json/include
                           external/pybind11_json/include
			   )

pybind11_add_module(gp_c gp.cpp )
#find_package(AMReX)
target_link_libraries(gp_c PRIVATE amrex ${MPI_LIBRARIES} gfortran gp_lib  )
#add_executable(main main.cpp run.cpp)
#target_link_libraries( main PUBLIC amrex ${MPI_LIBRARIES} gfortran pthread)

set (TEST_SOURCES tst/model-test.cpp)
add_executable(tst tst/main.cpp ${TEST_SOURCES})
target_link_libraries(tst PUBLIC  gtest  ${MPI_LIBRARIES} gp_lib amrex)
target_include_directories(tst PUBLIC
                           "${PROJECT_BINARY_DIR}/.."
                           )

install(TARGETS gp_c DESTINATION /home/luca/source/GP3/python)
