cmake_minimum_required(VERSION 3.10)
# set the project name
project(testBoundaries)
set(DIM 3)


set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_INSTALL_PREFIX /home/luca/source/GP3/install)

list(APPEND CMAKE_PREFIX_PATH CMAKE_PREFIX_PATH /home/luca/source/GP3/external/arrow_cpp)

add_subdirectory(../external/pybind11)
add_subdirectory(../external/amrex)
add_subdirectory(../external/json)
add_subdirectory(../external/pybind11_json)
add_subdirectory(../external/googletest)

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


set (TEST_SOURCES tst)
add_executable(tst tst/main.cpp ${TEST_SOURCES})
target_link_libraries(tst PRIVATE  gtest  ${MPI_LIBRARIES} amrex arrow parquet)
#target_link_libraries(tst PUBLIC  gp_lib)

target_include_directories(tst PUBLIC
                           "${PROJECT_BINARY_DIR}/.."
                           )