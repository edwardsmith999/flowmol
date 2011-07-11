#=======================================================================
# User definitions
#

NAME = setup

help:
	@echo The following targets are available:
	@echo
	@echo simple
	@echo blayer
	@echo channel
	@echo pipe
	@echo pipe.coax
	@echo jet.mask
	@echo jet.coax
	@echo jet.coax.dblexp
	@echo jet.coax.free
	@echo

SIMPLE  = main.o mesh.cartesian.o mesh.blayer.o initialField.blayer.o
BLAYER    = main.o mesh.cartesian.o mesh.blayer.o initialField.blayer.o RK4_Blasius.o Find_fpp_LengthInv.o Falkner_Skan.o
CHANNEL   = main.o mesh.cartesian.o mesh.channel.o initialField.channel.o
PIPE      = main.o mesh.cylindrical.o mesh.pipe.o initialField.pipe.o
PIPE_COAX = main.o mesh.cylindrical.o mesh.pipe.coax.o initialField.pipe.coax.o
JET_MASK  = main.o mesh.cylindrical.o mesh.jet.mask.o initialField.jet.coax.o
JET_COAX  = main.o mesh.cylindrical.o mesh.jet.coax.o initialField.jet.coax.o
JET_COAX_DBLEXP = main.o mesh.cylindrical.o mesh.jet.coax.dblexp.o initialField.jet.coax.o
JET_COAX_FREE   = main.o mesh.cylindrical.o mesh.jet.coax.free.o initialField.jet.coax.o

LIBFILES = $(HOME)/lib/archive.o $(HOME)/lib/linearMesh.o $(HOME)/lib/tanhMesh.o $(HOME)/lib/library.o $(HOME)/lib/library.alc.o

LIBS = 

FLAGS = $(MPI_COMPILE_FLAGS)

simple:
	@make opt "OFILES = $(SIMPLE)"
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
	@make all "FLAGS = -O3 $(FLAGS)"
clean::
	rm -rf obj a.out *.mod
a.out: obj $(OFILES)
	@echo "$(LD) $(OFILES) $(LIBFILES) $(LIBS) $(MPI_LD_FLAGS)"
	@cd obj; $(LD) $(OFILES) $(LIBFILES) $(LIBS) $(MPI_LD_FLAGS); mv a.out ..
obj:
	mkdir obj

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

