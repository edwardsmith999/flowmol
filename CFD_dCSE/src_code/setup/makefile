#=======================================================================
# User definitions
#
# Example of HECToR compilation line
# make LIB_BASE=./../../lib/hector-pgi USE_HITACHI=yes F90="ftn -r8" LD=ftn simple
#
#
#
NAME = setup

help:
	@echo The following targets are available:
	@echo
	@echo simple
	@echo read_in
	@echo blayer
	@echo channel
	@echo pipe
	@echo pipe.coax
	@echo jet.mask
	@echo jet.coax
	@echo jet.coax.dblexp
	@echo jet.coax.free
	@echo

SIMPLE    = mesh_parameters.o main.o mesh.cartesian.o mesh.blayer.o initialField.xz_periodic.o 
BLAYER    = mesh_parameters.o main.o mesh.cartesian.o mesh.blayer.o initialField.xz_periodic.o RK4_Blasius.o Find_fpp_LengthInv.o Falkner_Skan.o
CHANNEL   = mesh_parameters.o main.o mesh.cartesian.o mesh.channel.o initialField.channel.o
PIPE      = mesh_parameters.o main.o mesh.cylindrical.o mesh.pipe.o initialField.pipe.o
PIPE_COAX = mesh_parameters.o main.o mesh.cylindrical.o mesh.pipe.coax.o initialField.pipe.coax.o
JET_MASK  = mesh_parameters.o main.o mesh.cylindrical.o mesh.jet.mask.o initialField.jet.coax.o
JET_COAX  = mesh_parameters.o main.o mesh.cylindrical.o mesh.jet.coax.o initialField.jet.coax.o
JET_COAX_DBLEXP = mesh_parameters.o main.o mesh.cylindrical.o mesh.jet.coax.dblexp.o initialField.jet.coax.o
JET_COAX_FREE   = mesh_parameters.o main.o mesh.cylindrical.o mesh.jet.coax.free.o initialField.jet.coax.o

LIB_BASE := ./../../lib

LIBFILES := archive.o linearMesh.o tanhMesh.o library.o library.alc.o transforms.mkl_intel.alc.o
ifdef USE_HITACHI
  LIBFILES := archive.o linearMesh.o tanhMesh.o library.o library.alc.o transforms.fftw.hitachi.o
endif 

LIBFILES := $(addprefix $(LIB_BASE)/, $(LIBFILES))

#LIBFILES = $(HOME)/lib/archive.o $(HOME)/lib/linearMesh.o $(HOME)/lib/tanhMesh.o $(HOME)/lib/library.o $(HOME)/lib/library.alc.o $(HOME)/lib/SSDfile.memory.o $(HOME)/lib/triDiagonal.o $(HOME)/lib/triDiagonal.mpi.o $(HOME)/lib/pentaDiagonal.o $(HOME)/lib/transforms.mkl_intel.alc.o $(HOME)/lib/SinCosFFT.o $(HOME)/lib/triDiagonal.mrhs.o $(HOME)/lib/triDiagonal.mpi.mrhs.o $(HOME)/lib/pentaDiagonal.mrhs.o $(HOME)/lib/pentaDiagonal.mpi.o $(HOME)/lib/polyDiagonal.o

# LIBS = -L/apps/intel/cmkl/9.0/lib/64/ -lmkl_em64t -lmkl_dft_cluster -lguide -lpthread
# LIBS = -L/apps/intel/cmkl/9.0/lib/64/             -lmkl_dft_cluster -lguide -lpthread -lmkl
# LIBS = -L/apps/intel/cmkl/9.0/lib/64/ -lmkl_em64t                   -lguide -lpthread
# LIBS := -L/opt/intel/Compiler/11.1/046/mkl/lib/em64t/dns_lib -lmkl_em64t		  -lguide -lpthread
 LIBS := -mkl
#LIBS = -L/apps/intel/2013/mkl/lib/64/ -mkl

ifdef USE_HITACHI
  LIBS := -ldfftw -ldrfftw
endif

FLAGS = $(MPI_COMPILE_FLAGS)

simple:
	@make opt "OFILES = $(SIMPLE)"
read_in:
	@make opt "OFILES = $(READ_IN)"
blayer:
	@make opt "OFILES = $(BLAYER)"
channel:
	@make opt "OFILES = $(CHANNEL)"
pipe:
	@make opt "OFILES = $(PIPE)"
pipe.coax:
	@make opt "OFILES = $(PIPE_COAX)"
jet.mask:
	@make all "OFILES = $(JET_MASK)"
jet.coax:
	@make opt "OFILES = $(JET_COAX)"
jet.coax.dblexp:
	@make opt "OFILES = $(JET_COAX_DBLEXP)"
jet.coax.free:
	@make opt "OFILES = $(JET_COAX_FREE)"

#=======================================================================
# Standard definitions
#

# Cray compilers
# F90 = mpif90 -r8 -O3 -ip -funroll-loops -parallel -par-report3 -mcmodel=large -i-dynamic
F90 = mpif90 -r8 -O3 -ip -funroll-loops           -par-report3 -mcmodel=large -i-dynamic
#F90 = mpif90 -r8 -traceback -CB
F77 = f77
CC = cc
# LD = mpif90 -r8 -O3 -ip -funroll-loops -parallel -par-report3 -mcmodel=large -i-dynamic
LD = mpif90 -r8 -O3 -ip -funroll-loops           -par-report3 -mcmodel=large -i-dynamic

VPATH = obj
.SUFFIXES:
.SUFFIXES: .out .o .f90 .f .c .inc .h

#=======================================================================
# Targets and dependencies
#
default:: opt
all:: a.out
opt::
	@make all "FLAGS = $(FLAGS)"
clean::
	rm -rf obj a.out *.mod
a.out: obj $(OFILES)
	@echo "$(LD) $(OFILES) $(LIBFILES) $(LIBS) $(MPI_LD_FLAGS)"
	@cd obj; $(LD) $(OFILES) $(LIBFILES) $(LIBS) $(MPI_LD_FLAGS)  ; mv a.out ..
obj:
	mkdir -p obj

#=======================================================================
# Compilation rules
#
.f90.o:
	$(F90) $(FLAGS) -c $*.f90 -o obj/$*.o
.f.o:
	$(F77) $(FLAGS) -c $*.f
	@mv $*.o obj
.c.o:
	$(CC) $(FLAGS) -c $*.c -o obj/$*.o

