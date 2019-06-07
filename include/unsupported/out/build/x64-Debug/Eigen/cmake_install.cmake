# Install script for directory: E:/Users/Leoy/Documents/GitHub/Uovie-Hartree-Fock/lib/unsupported/Eigen

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "e:/users/leoy/documents/github/uovie-hartree-fock/lib/unsupported/out/install/x64-Debug")
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

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xDevelx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/unsupported/Eigen/AdolcForward;/unsupported/Eigen/AlignedVector3;/unsupported/Eigen/ArpackSupport;/unsupported/Eigen/AutoDiff;/unsupported/Eigen/BVH;/unsupported/Eigen/EulerAngles;/unsupported/Eigen/FFT;/unsupported/Eigen/IterativeSolvers;/unsupported/Eigen/KroneckerProduct;/unsupported/Eigen/LevenbergMarquardt;/unsupported/Eigen/MatrixFunctions;/unsupported/Eigen/MoreVectorization;/unsupported/Eigen/MPRealSupport;/unsupported/Eigen/NonLinearOptimization;/unsupported/Eigen/NumericalDiff;/unsupported/Eigen/OpenGLSupport;/unsupported/Eigen/Polynomials;/unsupported/Eigen/Skyline;/unsupported/Eigen/SparseExtra;/unsupported/Eigen/SpecialFunctions;/unsupported/Eigen/Splines")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/unsupported/Eigen" TYPE FILE FILES
    "E:/Users/Leoy/Documents/GitHub/Uovie-Hartree-Fock/lib/unsupported/Eigen/AdolcForward"
    "E:/Users/Leoy/Documents/GitHub/Uovie-Hartree-Fock/lib/unsupported/Eigen/AlignedVector3"
    "E:/Users/Leoy/Documents/GitHub/Uovie-Hartree-Fock/lib/unsupported/Eigen/ArpackSupport"
    "E:/Users/Leoy/Documents/GitHub/Uovie-Hartree-Fock/lib/unsupported/Eigen/AutoDiff"
    "E:/Users/Leoy/Documents/GitHub/Uovie-Hartree-Fock/lib/unsupported/Eigen/BVH"
    "E:/Users/Leoy/Documents/GitHub/Uovie-Hartree-Fock/lib/unsupported/Eigen/EulerAngles"
    "E:/Users/Leoy/Documents/GitHub/Uovie-Hartree-Fock/lib/unsupported/Eigen/FFT"
    "E:/Users/Leoy/Documents/GitHub/Uovie-Hartree-Fock/lib/unsupported/Eigen/IterativeSolvers"
    "E:/Users/Leoy/Documents/GitHub/Uovie-Hartree-Fock/lib/unsupported/Eigen/KroneckerProduct"
    "E:/Users/Leoy/Documents/GitHub/Uovie-Hartree-Fock/lib/unsupported/Eigen/LevenbergMarquardt"
    "E:/Users/Leoy/Documents/GitHub/Uovie-Hartree-Fock/lib/unsupported/Eigen/MatrixFunctions"
    "E:/Users/Leoy/Documents/GitHub/Uovie-Hartree-Fock/lib/unsupported/Eigen/MoreVectorization"
    "E:/Users/Leoy/Documents/GitHub/Uovie-Hartree-Fock/lib/unsupported/Eigen/MPRealSupport"
    "E:/Users/Leoy/Documents/GitHub/Uovie-Hartree-Fock/lib/unsupported/Eigen/NonLinearOptimization"
    "E:/Users/Leoy/Documents/GitHub/Uovie-Hartree-Fock/lib/unsupported/Eigen/NumericalDiff"
    "E:/Users/Leoy/Documents/GitHub/Uovie-Hartree-Fock/lib/unsupported/Eigen/OpenGLSupport"
    "E:/Users/Leoy/Documents/GitHub/Uovie-Hartree-Fock/lib/unsupported/Eigen/Polynomials"
    "E:/Users/Leoy/Documents/GitHub/Uovie-Hartree-Fock/lib/unsupported/Eigen/Skyline"
    "E:/Users/Leoy/Documents/GitHub/Uovie-Hartree-Fock/lib/unsupported/Eigen/SparseExtra"
    "E:/Users/Leoy/Documents/GitHub/Uovie-Hartree-Fock/lib/unsupported/Eigen/SpecialFunctions"
    "E:/Users/Leoy/Documents/GitHub/Uovie-Hartree-Fock/lib/unsupported/Eigen/Splines"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xDevelx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/unsupported/Eigen/src")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/unsupported/Eigen" TYPE DIRECTORY FILES "E:/Users/Leoy/Documents/GitHub/Uovie-Hartree-Fock/lib/unsupported/Eigen/src" FILES_MATCHING REGEX "/[^/]*\\.h$")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("E:/Users/Leoy/Documents/GitHub/Uovie-Hartree-Fock/lib/unsupported/out/build/x64-Debug/Eigen/CXX11/cmake_install.cmake")

endif()

