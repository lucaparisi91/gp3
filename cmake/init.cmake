set(DIM 3)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if ( ${CMAKE_BUILD_TYPE} MATCHES Debug)
  set(CMAKE_CXX_COMPILE_FLAGS_ADDITIONAL " -Wfatal-errors")
  set(CMAKE_CXX_LINK_FLAGS_ADDITIONAL "  ")

  
  elseif (${CMAKE_BUILD_TYPE} MATCHES Release) 
  set(CMAKE_CXX_COMPILE_FLAGS_ADDITIONAL " -DNDEBUG -Wfatal-errors -g -pg -march=core-avx2 -mfma")
  set(CMAKE_CXX_LINK_FLAGS_ADDITIONAL " ")

  else()
  message(FATAL_ERROR "Unrecognized build type: " ${CMAKE_BUILD_TYPE}  )

endif()

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

find_package(MPI REQUIRED)
find_package(OpenMP)
find_package(AMReX REQUIRED)

