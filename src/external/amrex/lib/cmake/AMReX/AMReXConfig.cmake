# ############################################################################ #
#
#  AMReX Configuration File
#  To import into other CMake projects
#
# ############################################################################ #

####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was AMReXConfig.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

# Set the minimum CMake version required
cmake_minimum_required(VERSION 3.14)

# Provides find_dependency
include(CMakeFindDependencyMacro)

#
# Build type
#
set(AMReX_BUILD_TYPE  Release)

#
# Versioning
#
set(AMReX_GIT_VERSION \"20.09-24-g82339abdc742-dirty\")

#
# AMReX CMake modules PATH
#
set_and_check(AMReX_MODULE_PATH ${PACKAGE_PREFIX_DIR}/Tools/CMake)

#
# Add AMReX modules to app code CMake
#
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${AMReX_MODULE_PATH})

#
# Configuration options
# Each option is treated like a "component" so that find_package can be easily
# used to check weather the option is enabled
#

# General options
set(AMReX_3D_FOUND                  ON)
set(AMReX_MPI_FOUND                 ON)
set(AMReX_MPI_THREAD_MULTIPLE_FOUND OFF)
set(AMReX_OMP_FOUND                 OFF)
set(AMReX_CUDA_FOUND                OFF)
set(AMReX_DPCPP_FOUND               OFF)
set(AMReX_DP_FOUND                  ON)

# Actual components selection
set(AMReX_EB_FOUND                  OFF)
set(AMReX_FINTERFACES_FOUND         OFF)
set(AMReX_LSOLVERS_FOUND            ON)
set(AMReX_AMRDATA_FOUND             OFF)
set(AMReX_PARTICLES_FOUND           OFF)
set(AMReX_DPARTICLES_FOUND          OFF)
set(AMReX_SENSEI_FOUND              OFF)
set(AMReX_CONDUIT_FOUND             OFF)
set(AMReX_SUNDIALS_FOUND            OFF)
set(AMReX_ASCENT_FOUND              OFF)
set(AMReX_HYPRE_FOUND               OFF)
set(AMReX_PETSC_FOUND               OFF)

# Compilation options
set(AMReX_FPE_FOUND                 OFF)
set(AMReX_PIC_FOUND                 OFF)
set(AMReX_ASSERTIONS_FOUND          OFF)

# Profiling options
set(AMReX_BASEP_FOUND               OFF)
set(AMReX_TINYP_FOUND               OFF)
set(AMReX_TRACEP_FOUND              OFF)
set(AMReX_MEMP_FOUND                OFF)
set(AMReX_COMMP_FOUND               OFF)
set(AMReX_PROFPARSER_FOUND          OFF)

#
# Parallel backends
#
set( THREADS_PREFER_PTHREAD_FLAG on)
find_dependency(Threads REQUIRED)

if (ON)
   set( _mpi_components C CXX )
   if (OFF)
      list(APPEND _mpi_components Fortran)
   endif ()
   find_dependency(MPI REQUIRED ${_mpi_components})
   unset(_mpi_components)
endif()

if (OFF)
   set( _omp_components CXX )
   if (OFF)
      list(APPEND _omp_components Fortran)
   endif ()
   find_dependency(OpenMP REQUIRED ${_omp_components})
endif ()

#
# Third party libraries
#
if (OFF)
   find_dependency(SENSEI REQUIRED)
endif ()

if (OFF)
    set(_sundials_components nvecserial;cvode;arkode)
    if (OFF)
        list(APPEND _sundials_components nvecopenmp)
    endif ()
    if (OFF)
        list(APPEND _sundials_components nveccuda)
    endif ()
    find_dependency(SUNDIALS 4 REQUIRED COMPONENTS ${_sundials_components})
    unset(_sundials_components)
endif ()

if (OFF)
    find_dependency(Ascent REQUIRED)
endif ()

if (OFF)
   find_dependency(Conduit REQUIRED)
endif ()

if (OFF)
   find_dependency(HYPRE 2.15 REQUIRED)
endif ()

if (OFF)
   find_dependency(PETSc 2.13 REQUIRED)
endif ()

#
# CUDA
#
if (OFF)
   include(AMReX_SetupCUDA)
endif ()


#include("${CMAKE_CURRENT_LIST_DIR}/.cmake")
include( "${CMAKE_CURRENT_LIST_DIR}/AMReXTargets.cmake" )

#
# Check components
#
check_required_components("AMReX")
#check_required_components(AMReX)
