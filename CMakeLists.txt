# -- This is the top-level CMakeLists.txt file for building semtex
#    codes and utilities.
#
# -- Set up is for out-of-source builds, typically run from within the
#    build subdirectory, e.g.: cd build; cmake ..; make

# -- Vendor-specific switches:

if (APPLE)

  set (BLA_VENDOR "Generic")
  message (STATUS "Building with Apple Accelerate BLAS and LAPACK.")

# -- Apple's C and C++ compilers may produce faster executables, but
#    the Gnu compilers gcc and g++ generate fewer warnings, link gnu's
#    gfortran libraries with little effort, and so are the default for
#    development on OS X.  To instead use XCode compilers, set the
#    following option to OFF.

  option (USE_GCC  "Build with gnu C & C++ compilers (not Xcode's)" OFF)

  if (USE_GCC) 
    set (CMAKE_C_COMPILER       "gcc")
    set (CMAKE_CXX_COMPILER     "g++")
    set (CMAKE_Fortran_COMPILER "gfortran")
  else ()
    # We have extra work to set RPATH and correctly link gfortran libraries.
    set (CMAKE_SKIP_BUILD_RPATH FALSE )
    if  (NOT DEFINED CMAKE_MACOSX_RPATH)
      set (CMAKE_MACOSX_RPATH ON)
    endif()
    # The following enables gcc12+fortran dylibs to be found on macports,
    # for now you should hand-edit the path if using homebrew or fink.
    set ( CMAKE_BUILD_RPATH "/opt/local/lib/libgcc" )
  endif (USE_GCC)
  
# -- N.B.: we have turned off all warnings via the -w flag.  Use with care!

  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 \
     -w ")
#    -Wno-parentheses \
#    -Wno-unused-value \
#    -Wno-format \
#    -Wno-format-security \
#    -Wno-switch \
#    -Wno-deprecated-declarations" )
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} \
     -w -std=c11 ")
#    -Wno-parentheses \
#    -Wno-implicit-function-declaration \
#    -Wno-format \
#    -Wno-format-security \
#    -Wno-deprecated-declarations" )
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=legacy" )
   
elseif (DEFINED ENV{BLA_VENDOR})

  # -- Typically this will indicate Intel compilers and MKL on a Linux system.
   
  message (STATUS "Building with vendor BLAS and LAPACK: " "$ENV{BLA_VENDOR}")

  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 \
    -w ")
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11 \
    -w ")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -f66" )

endif()

# -- End of vendor-specific setup.

cmake_minimum_required (VERSION 3.4)
project (semtex
	 VERSION     10.5
	 DESCRIPTION "Spectral element x Fourier DNS code"
	 LANGUAGES   C CXX Fortran )

# -- Pull in the pre-processor flag optons, which may also be set on
#    the comamnd line (if not enabled here).  E.g. cmake -DDEBUG=ON ..

option (USE_MPI  "Build dns+elliptic solvers with MPI"            ON)
option (BLD_DOG  "Build dog stability analysis solver and utils"  ON)
option (DEBUG    "Build with debugging preprocessor conditionals" OFF)

if (DEBUG)
  set (CMAKE_BUILD_TYPE "Debug")
  message (STATUS "Building code with debug preprocessing and compilation.")
else ()
  set (CMAKE_BUILD_TYPE "Release")
  message (STATUS "Building code with standard release optimizations.")
endif (DEBUG)

# -- External package dependencies.

find_package (BISON  REQUIRED)
find_package (BLAS   REQUIRED)
find_package (LAPACK REQUIRED)

if (USE_MPI)
  find_package (MPI REQUIRED)
  if (MPI_FOUND)
    message (STATUS "Building codes elliptic_mp and dns_mp with MPI support.")
  endif (MPI_FOUND)
endif (USE_MPI)

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

if (BLD_DOG)
  message (STATUS "Building codes for linear stability analysis.")
  include (dog/CMakeLists.txt)
endif (BLD_DOG)
