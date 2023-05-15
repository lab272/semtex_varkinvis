# -- Definitions used to compile elliptic solver.

set (ell_src
  ${CMAKE_CURRENT_LIST_DIR}/../src/message.cpp
  ${CMAKE_CURRENT_LIST_DIR}/helmholtz.cpp
  ${CMAKE_CURRENT_LIST_DIR}/drive.cpp
)

include_directories (veclib femlib ${MPI_C_INCLUDE_PATH})

set (import_serial_libs ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
set (import_par_libs    ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES}
    			${MPI_C_LIBRARIES})

# -- Default/serial executable (no MPI).

add_executable             (elliptic ${ell_src})
target_include_directories (elliptic PUBLIC src elliptic)
target_link_libraries      (elliptic src fem vec ${import_serial_libs})

# -- Add multi-process (_mp) executable if requested.

if (USE_MPI)
  add_executable             (elliptic_mp ${ell_src})
  target_compile_definitions (elliptic_mp PRIVATE -DMPI_EX)
  target_include_directories (elliptic_mp PUBLIC src elliptic)
  target_link_libraries      (elliptic_mp src fem vec ${import_par_libs})
endif (USE_MPI)

# -- Add partitioned (_xxt) executable if requested.

if (USE_XXT)
  add_executable             (elliptic_xxt ${ell_src})
  target_compile_definitions (elliptic_xxt PRIVATE -DXXT_EX)
  target_include_directories (elliptic_xxt PUBLIC src elliptic)
  target_link_libraries      (elliptic_xxt src fem vec ${import_par_libs})
endif (USE_XXT)
