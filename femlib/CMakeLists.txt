# -- NB:
#
#    The only semtex routine to directly require bison/yacc is initial.y.
#
#    The Fortran source files below are written in Fortran66/77.  Those
#    with extension .F will be run through the C preprocessor prior to
#    compilation; those with extension .f will be compiled without
#    preprocessing.

BISON_TARGET (Parser initial.y ${CMAKE_CURRENT_BINARY_DIR}/initial.c)

set (fem_lib_src
  ${BISON_Parser_OUTPUTS}
  polyops.c  operators.c   polylib.c    filter.c
  fourier.c  mapping.c     family.c
  temfftd.F  matops.F      sparsepak.F
  canfft.f   netlib.f  
)

add_library (fem STATIC ${fem_lib_src})
