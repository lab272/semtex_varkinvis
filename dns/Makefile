##############################################################################
# Makefile for in-source build of dns, spectral element solver for
# unsteady Navier--Stokes problems with generalised body forcing.
##############################################################################

PROG = dns

SEM   = ..
VPATH = $(SEM)/src

ifdef MPI		# -- Make a parallel executable.
  PROG := $(PROG)_mp
endif

# ----------------------------------------------------------------------------
# Semtex standard sources.
#
SEMFILES = analysis assemblymap auxfield bcmgr boundary boundarysys \
           condition domain  edge element family feml field geometry \
           history integration matrix mesh misc numbersys particle \
           statistics data2df svv
SEMOBJ   = $(addsuffix .o,$(SEMFILES))
SEMHDR   = $(addsuffix .h,$(SEMFILES)) sem.h

# ----------------------------------------------------------------------------
# Default build for Navier--Stokes solver.
#
NS_FILES = drive integrate integrates nonlinear dnsanalysis fieldforce
NS_OBJ   = $(addsuffix .o,$(NS_FILES))
NS_HDR   = dns.h fieldforce.h

$(PROG): sem $(SEMOBJ) $(NS_OBJ)
	$(LD) -o $@ $(NS_OBJ) $(SEMOBJ) $(LDFLAGS)

$(SEMOBJ): $(SEMHDR)
$(NS_OBJ): $(NS_HDR)

# ----------------------------------------------------------------------------
# Standard rules and definitions.
#
include $(VPATH)/Makefile

# ----------------------------------------------------------------------------
# Clean up executables and anything else
#
distclean:
	$(MAKE) clean
	rm -f dns dns_mp

#-----------------------------------------------------------------------------
# Machine-dependent compilations.

ifeq ($(ARCH),SUPER-UX)
mesh.o: $(VPATH)/mesh.cpp
	$(CXX) $(CPPFLAGS) $(NOVEC) $(VPATH)/mesh.cpp
endif
