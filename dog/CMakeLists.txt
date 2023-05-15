# if (EXISTS ${CMAKE_SOURCE_DIR}/dog)
#  option (WITH_FLIP   "Build dog with half-period map"                 OFF)
#  option (WITH_ARPACK "Build dog with ARPACK"                          OFF)
#  add_subdirectory (dog)

set (dog_src
     ${CMAKE_SOURCE_DIR}/dog/drive.cpp
     ${CMAKE_SOURCE_DIR}/dog/integrate.cpp
     ${CMAKE_SOURCE_DIR}/dog/stabanalysis.cpp

# -- Over-rides for standard semtex/src files:

     ${CMAKE_SOURCE_DIR}/dog/analysis.cpp
     ${CMAKE_SOURCE_DIR}/dog/auxfield.cpp
     ${CMAKE_SOURCE_DIR}/dog/bcmgr.cpp
     ${CMAKE_SOURCE_DIR}/dog/domain.cpp
     ${CMAKE_SOURCE_DIR}/dog/edge.cpp
     ${CMAKE_SOURCE_DIR}/dog/feml.cpp
     ${CMAKE_SOURCE_DIR}/dog/field.cpp
     ${CMAKE_SOURCE_DIR}/dog/geometry.cpp

# -- Standard semtex/src files:

     ${CMAKE_SOURCE_DIR}/src/assemblymap.cpp
     ${CMAKE_SOURCE_DIR}/src/boundary.cpp
     ${CMAKE_SOURCE_DIR}/src/boundarysys.cpp
     ${CMAKE_SOURCE_DIR}/src/condition.cpp
     ${CMAKE_SOURCE_DIR}/src/data2df.cpp		 
     ${CMAKE_SOURCE_DIR}/src/element.cpp
     ${CMAKE_SOURCE_DIR}/src/family.cpp
     ${CMAKE_SOURCE_DIR}/src/history.cpp
     ${CMAKE_SOURCE_DIR}/src/integration.cpp
     ${CMAKE_SOURCE_DIR}/src/matrix.cpp
     ${CMAKE_SOURCE_DIR}/src/mesh.cpp
     ${CMAKE_SOURCE_DIR}/src/misc.cpp
     ${CMAKE_SOURCE_DIR}/src/numbersys.cpp     
     ${CMAKE_SOURCE_DIR}/src/particle.cpp
     ${CMAKE_SOURCE_DIR}/src/svv.cpp
)

set (import_serial_libs ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

# -- Top-level executable (dog), default compilation.

add_executable             (dog ${dog_src})
target_include_directories (dog PUBLIC dog src)
target_link_libraries      (dog fem vec ${import_serial_libs})

# -- Optional executable (dog-AR), uses ARPACK for eigensystem solutions.

option (WITH_ARPACK "Build dog with ARPACK" OFF)
if (WITH_ARPACK)
  find_library (ARPACK_LIBRARY NAMES "arpack" HINTS ENV
  	       LD_LIBRARY_PATH PATHS
	       /usr/lib /usr/lib64 /usr/local/lib /usr/local/lib64)
  if (ARPACK_LIBRARY)
     message ("... found ARPACK: ${ARPACK_LIBRARY}")
  else()
     message (FATAL_ERROR "... have not found ARPACK, terminating build.")
  endif()
  add_executable             (dog-AR ${dog_src})
  target_compile_definitions (dog-AR PRIVATE ARPACK=ON)
  target_include_directories (dog-AR PUBLIC dog src)
  target_link_libraries      (dog-AR fem vec ${import_serial_libs}
  			     ${ARPACK_LIBRARY})	
endif (WITH_ARPACK)  

# option (WITH_FLIP   "Build dog with half-period map"                 OFF)
#  include (dog/BH-arnoldi/CMakeLists.txt)       
#  if (WITH_FLIP)
#    target_compile_definitions (dog PRIVATE FLIP=ON  )
#  endif (WITH_FLIP)
#  if (WITH_ARPACK)
#    target_compile_definitions (dog PRIVATE ARPACK=ON)
#  endif(WITH_ARPACK)

# -- lns is analogous to dns: it integrates the linearised Navier-Stokes
#    equations.  No eigensystem wrapper.
     
set (lns_src
     ${CMAKE_SOURCE_DIR}/dog/lns.cpp
     ${CMAKE_SOURCE_DIR}/dog/integrate.cpp
     ${CMAKE_SOURCE_DIR}/dog/stabanalysis.cpp

# -- Over-rides for standard semtex/src files:

     ${CMAKE_SOURCE_DIR}/dog/analysis.cpp
     ${CMAKE_SOURCE_DIR}/dog/auxfield.cpp
     ${CMAKE_SOURCE_DIR}/dog/bcmgr.cpp
     ${CMAKE_SOURCE_DIR}/dog/domain.cpp
     ${CMAKE_SOURCE_DIR}/dog/edge.cpp
     ${CMAKE_SOURCE_DIR}/dog/feml.cpp
     ${CMAKE_SOURCE_DIR}/dog/field.cpp
     ${CMAKE_SOURCE_DIR}/dog/geometry.cpp

# -- Standard semtex/src files:

     ${CMAKE_SOURCE_DIR}/src/assemblymap.cpp
     ${CMAKE_SOURCE_DIR}/src/boundary.cpp
     ${CMAKE_SOURCE_DIR}/src/boundarysys.cpp
     ${CMAKE_SOURCE_DIR}/src/condition.cpp
     ${CMAKE_SOURCE_DIR}/src/data2df.cpp		 
     ${CMAKE_SOURCE_DIR}/src/element.cpp
     ${CMAKE_SOURCE_DIR}/src/family.cpp
     ${CMAKE_SOURCE_DIR}/src/history.cpp
     ${CMAKE_SOURCE_DIR}/src/integration.cpp
     ${CMAKE_SOURCE_DIR}/src/matrix.cpp
     ${CMAKE_SOURCE_DIR}/src/mesh.cpp
     ${CMAKE_SOURCE_DIR}/src/misc.cpp
     ${CMAKE_SOURCE_DIR}/src/numbersys.cpp     
     ${CMAKE_SOURCE_DIR}/src/particle.cpp
     ${CMAKE_SOURCE_DIR}/src/svv.cpp
)

add_executable             (lns ${lns_src})
target_include_directories (lns PUBLIC dog src)
target_link_libraries      (lns fem vec ${import_serial_libs})

# -- Utilities follow.  See dog user guide or source file headers.
    
set (norm_src
     ${CMAKE_SOURCE_DIR}/dog/normalize.cpp
     ${CMAKE_SOURCE_DIR}/dog/auxfield.cpp
     ${CMAKE_SOURCE_DIR}/dog/feml.cpp
     ${CMAKE_SOURCE_DIR}/dog/geometry.cpp
     ${CMAKE_SOURCE_DIR}/src/data2df.cpp		 
     ${CMAKE_SOURCE_DIR}/src/element.cpp
     ${CMAKE_SOURCE_DIR}/src/family.cpp
     ${CMAKE_SOURCE_DIR}/src/mesh.cpp     
     ${CMAKE_SOURCE_DIR}/src/svv.cpp
)
add_executable             (normalize ${norm_src})
target_include_directories (normalize PUBLIC dog src)
target_link_libraries      (normalize fem vec ${import_serial_libs})

add_executable (symmetrise
     ${CMAKE_SOURCE_DIR}/dog/symmetrise.cpp
     ${CMAKE_SOURCE_DIR}/src/data2df.cpp
)     
target_include_directories (symmetrise PUBLIC dog src)
target_link_libraries      (symmetrise fem vec ${import_serial_libs})

add_executable (reflect
     ${CMAKE_SOURCE_DIR}/dog/reflect.cpp
     ${CMAKE_SOURCE_DIR}/src/data2df.cpp
)     
target_include_directories (reflect PUBLIC dog src)
target_link_libraries      (reflect fem vec ${import_serial_libs})

add_executable (combine ${CMAKE_SOURCE_DIR}/dog/combine.cpp)     
target_include_directories (combine PUBLIC dog)
target_link_libraries      (combine fem vec ${BLAS_LIBRARIES})

add_executable (flipmap ${CMAKE_SOURCE_DIR}/dog/flipmap.cpp)     
target_include_directories (flipmap PUBLIC dog src)
target_link_libraries      (flipmap fem vec ${BLAS_LIBRARIES})

configure_file (${CMAKE_SOURCE_DIR}/dog/circulate
	       	circulate.sh          COPYONLY)
