#############################################################################
# Makefile for lns, spectral element solver for linearised unsteady
# Navier--Stokes problems, and dog, a Krylov-method solver for
# stability problems (including Floquet and transient-growth problems).
#
##############################################################################

SEM   = ..
VPATH = $(SEM)/src

ifdef FLIP		# -- Use the "half-period-flip map".
  KPROG = dog-H
  DEFINES += -DFLIP
else
  KPROG = dog
endif

ifdef ARPACK
  DEFINES += -DARPACK
  KPROG   := $(KPROG)-AR
endif

# ----------------------------------------------------------------------------
# Semtex standard sources.
#
SEMFILES = analysis assemblymap auxfield bcmgr boundary boundarysys condition \
	   data2df domain edge element family feml field geometry history     \
           integration matrix mesh misc numbersys particle svv
SEMOBJ   = $(addsuffix .o,$(SEMFILES))
SEMHDR   = $(addsuffix .h,$(SEMFILES)) sem.h

# ----------------------------------------------------------------------------
# Build eigensolver code.

FL_FILES = drive integrate stabanalysis
FL_OBJ   = $(addsuffix .o,$(FL_FILES))

$(KPROG): sem $(FL_OBJ) $(SEMOBJ) stab.h
	$(LD) -o $@ $(FL_OBJ) $(SEMOBJ) $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build real-shift-invert eigensolver code.

SI_FILES = drive-rsi integrate stabanalysis templated
SI_OBJ   = $(addsuffix .o,$(SI_FILES))

dog-rsi: sem $(SI_OBJ) $(SEMOBJ) stab.h
	$(LD) -o $@ $(SI_OBJ) $(SEMOBJ) -L/usr/local/lib -larpack $(LDFLAGS)


# ----------------------------------------------------------------------------
# Build linearised Navier-Stokes solver.

ST_FILES = lns integrate stabanalysis
ST_OBJ   = $(addsuffix .o,$(ST_FILES))

lns: sem $(ST_OBJ) $(SEMOBJ) stab.h
	$(LD) -o $@ $(ST_OBJ) $(SEMOBJ) $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build combine, a utility to combine a base and a perturbation field.
#
combine: combine.o
	$(CXX) $(CPPFLAGS) -o $@ $@.cpp $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build normalise, a utility to normalise eigenvectors (2-norm based).
#
normalise: normalise.o
	$(CXX) $(CPPFLAGS) -o $@ $@.cpp $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build normalize, a utility to normalise eigenvectors (L2-norm based).
#

NRMSRC = normalize data2df element svv auxfield geometry family mesh feml
NRMOBJ   = $(addsuffix .o,$(NRMSRC))

normalize: $(NRMOBJ)
	$(CXX) $(CPPFLAGS) -o $@ $(NRMOBJ) $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build flipmap, a utility to find a symmetry-related mapping of the mesh.
#
flipmap: flipmap.cpp
	$(CXX) $(CPPFLAGS) -o $@ $@.cpp $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build symmetrise, a utility to enforce reflection symmstry of data.
#
symmetrise: symmetrise.cpp data2df.o 
	$(CXX) $(CPPFLAGS) -o $@ data2df.o $@.cpp $(LDFLAGS)

# ----------------------------------------------------------------------------
# Build reflect, a utility to generate a solution on a full mesh from
# one on a half mesh by reflection.
#
reflect: reflect.cpp data2df.o 
	$(CXX) $(CPPFLAGS) -o $@ data2df.o $@.cpp $(LDFLAGS)

# ----------------------------------------------------------------------------
# Header dependency.

$(SEMOBJ): $(SEMHDR)

# ----------------------------------------------------------------------------
# Standard rules and definitions.
#
include $(VPATH)/Makefile

ifdef ARPACK
  F77LIBS := -L/usr/local/lib -larpack $(F77LIBS)
endif

#-----------------------------------------------------------------------------
# Machine-dependent compilations.

ifeq ($(ARCH),SUPER-UX)
mesh.o: $(VPATH)/mesh.cpp
	$(CXX) $(CPPFLAGS) $(NOVEC) $(VPATH)/mesh.cpp
endif
