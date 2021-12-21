# -- This is the top-level CMakeLists.txt file for building semtex
#    codes and utilities.
#
# -- Set up is for out-of-source builds, typically run from within the
#    build subdirectory, e.g.: cd build; cmake ..; make

project (semtex C CXX Fortran)
cmake_minimum_required (VERSION 3.0)

# -- Static linkage (comment out for dynamic linkage/default).

#SET (CMAKE_FIND_LIBRARY_SUFFIXES ".a")
#SET (BUILD_SHARED_LIBS OFF)
#SET (CMAKE_EXE_LINKER_FLAGS "-static")

# -- Pull in the pre-processor flags, which can optionally be set on
#    the comamnd line (if not enabled here).  E.g. cmake -DDEBUG=ON ..
#
# -- Note that even with the MPI option set, all codes other than the
#    elliptic and dns solvers are serial executables and hence do not
#    exploit or benefit from presence of MPI.

option (WITH_MPI    "Build dns+elliptic solvers with MPI"            ON )
option (DEBUG       "Build with debugging preprocessor conditionals" OFF)
if (DEBUG)
  set (CMAKE_BUILD_TYPE "Debug")
  message (STATUS "Building code with debug preprocessing and compilation.")
else ()
  set (CMAKE_BUILD_TYPE "Release")
  message (STATUS "Building code with standard release optimizations.")
endif ()

# -- Vendor-specific switches:

if (${CMAKE_HOST_SYSTEM_NAME} STREQUAL "Darwin")
   set (BLA_VENDOR "Generic")
   message (STATUS "Building with Apple Accelerate BLAS and LAPACK.")
   set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} \
     -Wno-parentheses \
     -Wno-unused-value \
     -Wno-implicit-function-declaration \
     -Wno-format \
     -Wno-format-security -Wno-switch")
   set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} \
     -Wno-parentheses \
     -Wno-format \
     -Wno-format-security")
elseif (DEFINED ENV{BLA_VENDOR})
   message (STATUS "Building with vendor BLAS and LAPACK: " "$ENV{BLA_VENDOR}")
endif()	   

# -- External package dependencies.

find_package   (BISON  REQUIRED)
find_package   (BLAS   REQUIRED)
find_package   (LAPACK REQUIRED)
if (WITH_MPI)
  find_package (MPI)
  if (MPI_FOUND)
    message (STATUS "Building codes elliptic_mp and dns_mp with MPI support.")
  endif (MPI_FOUND)
endif (WITH_MPI)

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
