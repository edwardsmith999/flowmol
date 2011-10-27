#=======================================================================
# User definitions
#

NAME = library
LIBDIR = ../lib

F90FILES = archive.f90 library.f90 mkl_dfti.f90 pentaDiagonal.f90  SinCosFFT.f90  timer.f90  triDiagonal.mpi.mrhs.f90  triDiagonal.f90  linearMesh.f90  pentaDiagonal.mpi.f90  polyDiagonal.f90  SSDfile.memory.f90  timeStep.RK3.f90  triDiagonal.mpi.f90  library.alc.f90  math.f90  pentaDiagonal.mrhs.f90  quadrature.f90  tanhMesh.f90  transforms.mkl_intel.alc.f90 triDiagonal.mrhs.f90

#  LIBS     = #-L$(MKL_HOME)                 -lmkl_em64t
 LIBS     = -L/apps/intel/cmkl/9.0/lib/64 -lmkl_em64t

FLAGS =
F77FILES = 
CFILES = 

#=======================================================================
# Standard definitions
#
OFILES = $(F90FILES:.f90=.o) $(F77FILES:.f=.o) $(CFILES:.c=.o)

COMP := GNU

# INTEL compilers
# F90 = mpif90 -r8 -O3 -ip -funroll-loops -parallel -mcmodel=large -i-dynamic
  F90_INTEL = mpif90 -r8 -O3 -ip -funroll-loops           -mcmodel=large -i-dynamic
  F77_INTEL = f77
  CC_INTEL = cc
# LD = mpif90  -r8 -O3 -ip -funroll-loops -parallel -mcmodel=large -i-dynamic
  LD_INTEL := mpif90  -r8 -O3 -ip -funroll-loops           -mcmodel=large -i-dynamic
 
# HECToR GNU compiler

  F90_GNU := ftn -O3 -fdefault-real-8 -fdefault-double-8
  F77_GNU := ftn
  CC_GNU  := cc
  LD_GNU  := ftn -O3 -fdefault-real-8 -fdefault-double-8

# HECToR PGI compiler

  F90_PGI := ftn -fast -r8
  F77_PGI := ftn -fast
  CC_PGI  := cc
  LD+PGI  := ftn -r8

F90 := $(F90_$(COMP))
F77 := $(F77_$(COMP))
CC  := $(cc_$(COMP))
LD  := $(LD_$(COMP))


VPATH = obj
.SUFFIXES:
.SUFFIXES: .out .o .f90 .f .c .inc .h

#=======================================================================
# Targets and dependencies
#
default:: opt
all::   $(LIBDIR) $(OFILES)
	rm -f *.mod
clean:
	rm -f *.mod
	cd $(LIBDIR)
	rm -f $(OFILES)
debug:
	@make all "FLAGS = -g $(FLAGS)"
prof:
	@make opt "LIBS = -p $(LIBS)"
opt:
	@make all "FLAGS = -O3 $(FLAGS)"
$(LIBDIR):
	[ -d $(LIBDIR) ] || mkdir $(LIBDIR)

#=======================================================================
# Compilation rules
#
.f90.o:
	$(F90) $(FLAGS) -c $*.f90 -o $(LIBDIR)/$*.o
.f.o:
	$(F77) $(FLAGS) -c $*.f -o $*.o
.c.o:
	$(CC)  $(FLAGS) -c $*.c -o $*.o

