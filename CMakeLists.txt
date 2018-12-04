# -- This is the top-level CMakeLists.txt file for building semtex
#    codes and utilities.
#
# -- Set up is for out-of-source builds, typically run from within the
#    build subdirectory, e.g.: cd build; cmake ..; make

project (semtex C CXX Fortran)
cmake_minimum_required (VERSION 2.6)

# -- Pull in the pre-processor flags, which can optionally be set on
#    the comamnd line (if not enabled here).  E.g. cmake -DWITH_MPI=ON ..
#
# -- Note that even with the MPI option set, all codes other than the
#    elliptic and dns solvers are serial executables and hence do not
#    exploit or benefit from presence of MPI.

option (WITH_MPI    "Build dns+elliptic solvers with MPI"            OFF)
option (DEBUG       "Build with debugging preprocessor conditionals" OFF)
if (DEBUG)
  set (CMAKE_BUILD_TYPE "Debug")
  message (STATUS "Building code with debug preprocessing and compilation.")
else ()
  set (CMAKE_BUILD_TYPE "Release")
  message (STATUS "Building code with standard release optimizations.")
endif ()

# -- External package dependencies.

find_package   (BISON  REQUIRED)
find_package   (BLAS   REQUIRED)
find_package   (LAPACK REQUIRED)
if (WITH_MPI)
  find_package (MPI    REQUIRED)
endif (WITH_MPI)

# -- Set compiler warning supressions (some of which are unknown to Fortran).

add_compile_options (
  -Wno-parentheses
  -Wno-unused-value
  -Wno-implicit-function-declaration
  -Wno-format
  -Wno-format-security  
  -Wno-switch
)

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

if (EXISTS ${CMAKE_SOURCE_DIR}/stab)
  include (${CMAKE_SOURCE_DIR}/stab/CMakeLists.txt)
endif()
