

set(PYTHON_BINDING_FILES gpIO.cpp saveMultifab.cpp )
set( PYBIND11_FINDPYTHON ON )

find_package(PythonLibs)


include_directories(
${PROJECT_SOURCE_DIR}/external/json/include
${PROJECT_SOURCE_DIR}/external/pybind11_json/include
			   )

include_directories(
${PROJECT_SOURCE_DIR}/external/json/include
			   )


pybind11_add_module(gpIO_c ${PYTHON_BINDING_FILES} )
include_directories( ${PYTHON_INCLUDE_DIRS})
target_link_libraries( gpIO_c PRIVATE amrex ${MPI_LIBRARIES} gfortran gp_lib )

