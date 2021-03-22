# Install script for directory: /home/luca/source/GP3/1D/external/eigen/Eigen

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xDevelx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/Eigen" TYPE FILE FILES
    "/home/luca/source/GP3/1D/external/eigen/Eigen/Cholesky"
    "/home/luca/source/GP3/1D/external/eigen/Eigen/CholmodSupport"
    "/home/luca/source/GP3/1D/external/eigen/Eigen/Core"
    "/home/luca/source/GP3/1D/external/eigen/Eigen/Dense"
    "/home/luca/source/GP3/1D/external/eigen/Eigen/Eigen"
    "/home/luca/source/GP3/1D/external/eigen/Eigen/Eigenvalues"
    "/home/luca/source/GP3/1D/external/eigen/Eigen/Geometry"
    "/home/luca/source/GP3/1D/external/eigen/Eigen/Householder"
    "/home/luca/source/GP3/1D/external/eigen/Eigen/IterativeLinearSolvers"
    "/home/luca/source/GP3/1D/external/eigen/Eigen/Jacobi"
    "/home/luca/source/GP3/1D/external/eigen/Eigen/LU"
    "/home/luca/source/GP3/1D/external/eigen/Eigen/MetisSupport"
    "/home/luca/source/GP3/1D/external/eigen/Eigen/OrderingMethods"
    "/home/luca/source/GP3/1D/external/eigen/Eigen/PaStiXSupport"
    "/home/luca/source/GP3/1D/external/eigen/Eigen/PardisoSupport"
    "/home/luca/source/GP3/1D/external/eigen/Eigen/QR"
    "/home/luca/source/GP3/1D/external/eigen/Eigen/QtAlignedMalloc"
    "/home/luca/source/GP3/1D/external/eigen/Eigen/SPQRSupport"
    "/home/luca/source/GP3/1D/external/eigen/Eigen/SVD"
    "/home/luca/source/GP3/1D/external/eigen/Eigen/Sparse"
    "/home/luca/source/GP3/1D/external/eigen/Eigen/SparseCholesky"
    "/home/luca/source/GP3/1D/external/eigen/Eigen/SparseCore"
    "/home/luca/source/GP3/1D/external/eigen/Eigen/SparseLU"
    "/home/luca/source/GP3/1D/external/eigen/Eigen/SparseQR"
    "/home/luca/source/GP3/1D/external/eigen/Eigen/StdDeque"
    "/home/luca/source/GP3/1D/external/eigen/Eigen/StdList"
    "/home/luca/source/GP3/1D/external/eigen/Eigen/StdVector"
    "/home/luca/source/GP3/1D/external/eigen/Eigen/SuperLUSupport"
    "/home/luca/source/GP3/1D/external/eigen/Eigen/UmfPackSupport"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xDevelx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/Eigen" TYPE DIRECTORY FILES "/home/luca/source/GP3/1D/external/eigen/Eigen/src" FILES_MATCHING REGEX "/[^/]*\\.h$")
endif()

