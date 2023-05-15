# -- This is the top-level CMakeLists.txt file for building semtex
#    codes and utilities.
#
# -- Set up is for out-of-source builds, typically run from within the
#    build subdirectory, e.g.: cd build; cmake ..; make

cmake_minimum_required (VERSION 3.0)

# -- Vendor-specific switches:

if (${CMAKE_HOST_SYSTEM_NAME} STREQUAL "Darwin")

# -- Using Apple's C and C++ compilers may produce faster executables,
#    but the choices gcc and g++ below generate fewer warnings, and so
#    are the default for development on OS X.  Comment out these
#    compiler definitions for released versions:

#  set (CMAKE_C_COMPILER       "gcc")
#  set (CMAKE_CXX_COMPILER     "g++")
#  set (CMAKE_Fortran_COMPILER "gfortran")

  set (BLA_VENDOR "Generic")
  message (STATUS "Building with Apple Accelerate BLAS and LAPACK.")

# -- N.B.: we have turned off all warnings via the -w flag.  Use with care!

  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} \
    -w \
    -Wno-parentheses \
    -Wno-unused-value \
    -Wno-format \
    -Wno-format-security \
    -Wno-switch \
    -Wno-deprecated-declarations" )
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} \
    -w \
    -Wno-parentheses \
    -Wno-implicit-function-declaration \
    -Wno-format \
    -Wno-format-security \
    -Wno-deprecated-declarations" )
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=legacy" )
   
elseif (DEFINED ENV{BLA_VENDOR})

  # -- Typically this will indicate Intel compilers and MKL.
   
  message (STATUS "Building with vendor BLAS and LAPACK: " "$ENV{BLA_VENDOR}")

  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} \
    -w \
    -std=c++11 \
    -Wno-parentheses \
    -Wno-unused-value \
    -Wno-format \
    -Wno-format-security \
    -Wno-switch \
    -Wno-deprecated-declarations" )
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} \
    -w \
    -std=c11 \
    -Wno-parentheses \
    -Wno-implicit-function-declaration \
    -Wno-format \
    -Wno-format-security \
    -Wno-deprecated-declarations" )
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -f66" )

endif()	   

project (semtex C CXX Fortran)

# -- Pull in the pre-processor flags, which can optionally be set on
#    the comamnd line (if not enabled here).  E.g. cmake -DDEBUG=ON ..

option (USE_MPI  "Build dns+elliptic solvers with MPI"            ON )
option (USE_XXT  "Build dns+elliptic solvers with XXT"            ON )
option (DEBUG    "Build with debugging preprocessor conditionals" OFF)

if (DEBUG)
  set (CMAKE_BUILD_TYPE "Debug")
  message (STATUS "Building code with debug preprocessing and compilation.")
else ()
  set (CMAKE_BUILD_TYPE "Release")
  message (STATUS "Building code with standard release optimizations.")
endif (DEBUG)

# -- External package dependencies.

find_package   (BISON  REQUIRED)
find_package   (BLAS   REQUIRED)
find_package   (LAPACK REQUIRED)

if (USE_MPI)
  find_package (MPI REQUIRED)
  if (MPI_FOUND)
    message (STATUS "Building codes elliptic_mp and dns_mp with MPI support.")
  endif (MPI_FOUND)
endif (USE_MPI)

# -- For 2D partitioned solver option, we will use XXT and Scotch together.

if (USE_XXT)

  if (NOT MPI_FOUND)
    message (FATAL_ERROR "MPI is needed to compile with XXT.")
  endif (NOT MPI_FOUND)

  # -- Set up to compile third-party library dependencies that we can
  #    supply, as they are less likely to be already available.

  #    Add support for CMAKE_DEPENDENT_OPTION:

  include (CMakeDependentOption)
  include (CMakeParseArguments)

  set (CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})
  set (TPSRC   ${CMAKE_SOURCE_DIR}/ThirdParty)
  set (TPBUILD ${CMAKE_SOURCE_DIR}/ThirdParty)
  set (TPDIST  ${CMAKE_BINARY_DIR}/ThirdParty/dist)
  set (TPINC   ${CMAKE_INSTALL_PREFIX}/include/ThirdParty)
  # -- the following for monash storage:
  #set (TPURL   https://users.monash.edu.au/~bburn/src)
  # -- the following for gitlab storage:
  # -- see https://stackoverflow.com/questions/56943327
  set (TPURL https://gitlab.com/api/v4/projects/32564913/repository/files)

  include (SemtexCommon)
  include (ThirdPartyScotch)
  include (ThirdPartyMPI)

endif (USE_XXT)

# -- Build semtex libraries.

include_directories (veclib femlib ${MPI_C_INCLUDE_PATH})

add_subdirectory (veclib)
add_subdirectory (femlib)
add_subdirectory (src)

# -- Build the top-level spectral element solvers elliptic and dns.

include (elliptic/CMakeLists.txt)
include (dns/CMakeLists.txt)

# -- Build/install pre/post processing utilities.

include (utility/CMakeLists.txt)

# -- Regression tests for dns and elliptic and their parallel
#    equivalents (the latter, if requested).  After building
#    completes, please run ctest manually within the build directory.

enable_testing ()
include (test/CMakeLists.txt)

# -- Optional compilation of linear stability codes, if supplied.

if (EXISTS ${CMAKE_SOURCE_DIR}/dog)
  include (${CMAKE_SOURCE_DIR}/dog/CMakeLists.txt)
endif()
