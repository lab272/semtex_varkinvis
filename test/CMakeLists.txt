# -- Serial tests of elliptic (the same elliptic solver used in dns).

add_test(laplace1 ${CMAKE_SOURCE_DIR}/test/testregression ""
		 ${CMAKE_SOURCE_DIR}/build/elliptic laplace1)
add_test(laplace3 ${CMAKE_SOURCE_DIR}/test/testregression ""
		 ${CMAKE_SOURCE_DIR}/build/elliptic laplace3)
add_test(laplace4 ${CMAKE_SOURCE_DIR}/test/testregression ""
		 ${CMAKE_SOURCE_DIR}/build/elliptic laplace4)
add_test(laplace5 ${CMAKE_SOURCE_DIR}/test/testregression ""
		 ${CMAKE_SOURCE_DIR}/build/elliptic laplace5)
add_test(laplace6 ${CMAKE_SOURCE_DIR}/test/testregression ""
		 ${CMAKE_SOURCE_DIR}/build/elliptic laplace6)
add_test(laplace7 ${CMAKE_SOURCE_DIR}/test/testregression ""
		 ${CMAKE_SOURCE_DIR}/build/elliptic laplace7)
add_test(poisson1 ${CMAKE_SOURCE_DIR}/test/testregression ""
		 ${CMAKE_SOURCE_DIR}/build/elliptic poisson1)
add_test(poisson2 ${CMAKE_SOURCE_DIR}/test/testregression ""
		 ${CMAKE_SOURCE_DIR}/build/elliptic poisson2)
add_test(helmholtz1 ${CMAKE_SOURCE_DIR}/test/testregression ""
		 ${CMAKE_SOURCE_DIR}/build/elliptic helmholtz1)
add_test(mixed1 ${CMAKE_SOURCE_DIR}/test/testregression ""
		 ${CMAKE_SOURCE_DIR}/build/elliptic mixed1)

# -- Serial tests of dns for both 2D and 3D problems:

add_test(taylor2 ${CMAKE_SOURCE_DIR}/test/testregression ""
		 ${CMAKE_SOURCE_DIR}/build/dns taylor2)
add_test(taylor3 ${CMAKE_SOURCE_DIR}/test/testregression ""
		 ${CMAKE_SOURCE_DIR}/build/dns taylor3)
add_test(taylor4 ${CMAKE_SOURCE_DIR}/test/testregression ""
		 ${CMAKE_SOURCE_DIR}/build/dns taylor4)
add_test(taylor5 ${CMAKE_SOURCE_DIR}/test/testregression ""
		 ${CMAKE_SOURCE_DIR}/build/dns taylor5)
add_test(kovas1  ${CMAKE_SOURCE_DIR}/test/testregression ""
		 ${CMAKE_SOURCE_DIR}/build/dns kovas1 )
add_test(kovas2  ${CMAKE_SOURCE_DIR}/test/testregression ""
		 ${CMAKE_SOURCE_DIR}/build/dns kovas2 )
add_test(kovas3  ${CMAKE_SOURCE_DIR}/test/testregression ""
		 ${CMAKE_SOURCE_DIR}/build/dns kovas3 )
add_test(kovas4  ${CMAKE_SOURCE_DIR}/test/testregression ""
		 ${CMAKE_SOURCE_DIR}/build/dns kovas4 )
add_test(kovas5  ${CMAKE_SOURCE_DIR}/test/testregression ""
		 ${CMAKE_SOURCE_DIR}/build/dns kovas5 )
add_test(tube1   ${CMAKE_SOURCE_DIR}/test/testregression ""
		 ${CMAKE_SOURCE_DIR}/build/dns tube1  )
add_test(tube2   ${CMAKE_SOURCE_DIR}/test/testregression ""
		 ${CMAKE_SOURCE_DIR}/build/dns tube2  )
add_test(tube3   ${CMAKE_SOURCE_DIR}/test/testregression ""
		 ${CMAKE_SOURCE_DIR}/build/dns tube3  )
add_test(sbr     ${CMAKE_SOURCE_DIR}/test/testregression ""
		 ${CMAKE_SOURCE_DIR}/build/dns sbr    )
add_test(tc1     ${CMAKE_SOURCE_DIR}/test/testregression ""
		 ${CMAKE_SOURCE_DIR}/build/dns tc1    )
add_test(cylkov2 ${CMAKE_SOURCE_DIR}/test/testregression ""
		 ${CMAKE_SOURCE_DIR}/build/dns cylkov2)

# -- Parallel tests of elliptic and dns for 3D problems:

if(WITH_MPI)
  add_test(laplace3_mp ${CMAKE_SOURCE_DIR}/test/testregression "mpirun -np 2"
  		   ${CMAKE_SOURCE_DIR}/build/elliptic_mp laplace3)
  add_test(taylor3_mp ${CMAKE_SOURCE_DIR}/test/testregression "mpirun -np 2"
  		   ${CMAKE_SOURCE_DIR}/build/dns_mp taylor3)
  add_test(taylor4_mp ${CMAKE_SOURCE_DIR}/test/testregression "mpirun -np 2"
  		   ${CMAKE_SOURCE_DIR}/build/dns_mp taylor4)
  add_test(taylor5_mp ${CMAKE_SOURCE_DIR}/test/testregression "mpirun -np 2"
  		   ${CMAKE_SOURCE_DIR}/build/dns_mp taylor5)
  add_test(kovas2_mp  ${CMAKE_SOURCE_DIR}/test/testregression "mpirun -np 2"
  		   ${CMAKE_SOURCE_DIR}/build/dns_mp kovas2 )
  add_test(kovas3_mp  ${CMAKE_SOURCE_DIR}/test/testregression "mpirun -np 2"
  		   ${CMAKE_SOURCE_DIR}/build/dns_mp kovas3 )
  add_test(kovas4_mp  ${CMAKE_SOURCE_DIR}/test/testregression "mpirun -np 2"
  		   ${CMAKE_SOURCE_DIR}/build/dns_mp kovas4 )
  add_test(kovas5_mp  ${CMAKE_SOURCE_DIR}/test/testregression "mpirun -np 2"
  		   ${CMAKE_SOURCE_DIR}/build/dns_mp kovas5 )
endif()
