# -- Definitions used to compile linear stability analysis solver
#    "dog" and associated utilities.

# -- Refer to drive.cpp regarding these options for dog:

option (WITH_FLIP   "Build dog with half-period map"  OFF )
option (WITH_ARPACK "Build dog with ARPACK"           Off )

# -- Local overrides for standard semtex/src files:

set (override_src
     ${CMAKE_CURRENT_LIST_DIR}/analysis.cpp
     ${CMAKE_CURRENT_LIST_DIR}/auxfield.cpp
     ${CMAKE_CURRENT_LIST_DIR}/bcmgr.cpp
     ${CMAKE_CURRENT_LIST_DIR}/domain.cpp
     ${CMAKE_CURRENT_LIST_DIR}/edge.cpp
     ${CMAKE_CURRENT_LIST_DIR}/feml.cpp
     ${CMAKE_CURRENT_LIST_DIR}/field.cpp
     ${CMAKE_CURRENT_LIST_DIR}/geometry.cpp
)
add_library (override STATIC ${override_src})
target_include_directories (override PUBLIC dog src)

# -- Top-level executable (dog), default compilation.

set (dog_src
     ${CMAKE_CURRENT_LIST_DIR}/drive.cpp
     ${CMAKE_CURRENT_LIST_DIR}/integrate.cpp
     ${CMAKE_CURRENT_LIST_DIR}/stabanalysis.cpp
)
add_executable             (dog ${dog_src})
target_include_directories (dog PUBLIC dog src)
target_link_libraries      (dog override src fem vec ${import_serial_libs})

# -- Optional executable dog-AR, uses ARPACK for eigensystem solutions.

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
  target_link_libraries      (dog-AR override src fem vec ${import_serial_libs}
  			     ${ARPACK_LIBRARY})	

endif(WITH_ARPACK)

# -- Optional executable dog-H, incorporates half-period-flip map.

if (WITH_FLIP)
  add_executable             (dog-H ${dog_src})
  target_compile_definitions (dog-H PRIVATE FLIP=ON)
  target_include_directories (dog-H PUBLIC dog src)
  target_link_libraries      (dog-H override src fem vec ${import_serial_libs})
endif (WITH_FLIP)

# -- lns is analogous to dns: it integrates the linearised Navier-Stokes
#    equations.  No eigensystem wrapper.
     
set (lns_src
     ${CMAKE_CURRENT_LIST_DIR}/lns.cpp
     ${CMAKE_CURRENT_LIST_DIR}/integrate.cpp
     ${CMAKE_CURRENT_LIST_DIR}/stabanalysis.cpp
)
add_executable             (lns ${lns_src})
target_include_directories (lns PUBLIC dog src)
target_link_libraries      (lns override src fem vec ${import_serial_libs})

# -- Utilities follow.  See dog user guide or source file headers.
    
set (norm_src
     ${CMAKE_CURRENT_LIST_DIR}/normalize.cpp
)
add_executable            (normalize ${norm_src})
target_include_directories(normalize PUBLIC dog src)
target_link_libraries     (normalize override src fem vec ${import_serial_libs})

add_executable (symmetrise
     ${CMAKE_CURRENT_LIST_DIR}/symmetrise.cpp
     ${CMAKE_CURRENT_LIST_DIR}/../src/data2df.cpp
)     
target_include_directories (symmetrise PUBLIC dog src)
target_link_libraries      (symmetrise fem vec ${import_serial_libs})

add_executable (reflect
     ${CMAKE_CURRENT_LIST_DIR}/reflect.cpp
     ${CMAKE_CURRENT_LIST_DIR}/../src/data2df.cpp
)     
target_include_directories (reflect PUBLIC dog src)
target_link_libraries      (reflect fem vec ${import_serial_libs})

add_executable (combine ${CMAKE_CURRENT_LIST_DIR}/combine.cpp)     
target_include_directories (combine PUBLIC dog)
target_link_libraries      (combine fem vec ${BLAS_LIBRARIES})

add_executable (flipmap ${CMAKE_CURRENT_LIST_DIR}/flipmap.cpp)     
target_include_directories (flipmap PUBLIC dog src)
target_link_libraries      (flipmap fem vec ${BLAS_LIBRARIES})

configure_file (${CMAKE_CURRENT_LIST_DIR}/circulate circulate.sh COPYONLY)
