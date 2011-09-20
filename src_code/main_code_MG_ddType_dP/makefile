#=======================================================================
# User definitions
#

NAME = blayer

F90FILES := main.f90 cfd_control.f90 data.export_blocks.f90 mesh.export.f90 grid.export.f90 data.noglobal_blocks.f90 arrayzero.f90 messenger.MPI.noglobal_blocks.f90 misc.mpi.f90 mesh.f90 subMesh.f90 grid.f90 mtgcoeff_blocks.f90 inflow.f90 outflow.f90 globalBCs.f90 fluct_random.f90 fluctuations_inlet_top.f90 advance_quick_Fletcher_xQUICK_deltap.f90 invert_sgi_openmp.f90 poisson_newmtg_blocks.sgi_deltap.f90 multigrid_blocks_ylines_zebra_FullMG.f90 projection_step_deltap.f90 cfl_control.f90 copyBC.f90 flowAnalysis.f90 statistics.f90 Output_SubDomain.f90 Input_mod.Z_blocks.f90 Output_mod.Z_blocks.f90 Input_Spectra.f90 Output_Spectra.f90 Input_Stats.f90 Output_Stats.f90 simulation.f90 coupler_cfd_setup.f90 coupler_cfd_global_data.f90 coupler_cfd_communication.f90

# use coupler ? 
USE_COUPLER      := no
COUPLER_PATH_yes := ../../../coupler_dCSE/src_code/
COUPLER_PATH_no  := ./
COUPLER_PATH     := $(COUPLER_PATH_$(USE_COUPLER))


F77FILES = 
CFILES = 

LIBBASE=../lib

LIBFILES = $(LIBBASE)/archive.o $(LIBBASE)/library.o $(LIBBASE)/library.alc.o $(LIBBASE)/SSDfile.memory.o $(LIBBASE)/triDiagonal.o $(LIBBASE)/triDiagonal.mpi.o $(LIBBASE)/pentaDiagonal.o $(LIBBASE)/transforms.mkl_intel.alc.o $(LIBBASE)/SinCosFFT.o $(LIBBASE)/triDiagonal.mrhs.o $(LIBBASE)/triDiagonal.mpi.mrhs.o $(LIBBASE)/pentaDiagonal.mrhs.o $(LIBBASE)/pentaDiagonal.mpi.o $(LIBBASE)/polyDiagonal.o 

#  LIBS = -L/apps/intel/cmkl/9.0/lib/64/ -lmkl_em64t -lmkl_dft_cluster -lguide -lpthread
#  LIBS = -L/opt/intel/cmkl/9.0/lib/64/             -lmkl_dft_cluster -lguide -lpthread -lmkl
#  LIBS = -L/apps/intel/cmkl/9.0/lib/64/ -lmkl_em64t                   -lguide -lpthread
LIBS = -L/opt/intel/Compiler/11.1/046/mkl/lib/em64t/ -lmkl_em64t                   -lguide -lpthread


#=======================================================================
# Compiler
#=======================================================================

COMP    = gnu
F90     = mpif90
F77     = mpif90
CC      = cc
LD      = mpif90
VERSION = opt

# GNU flags
FLAGS_gnu_debug = -fdefault-real-8 -fdefault-double-8 -ffree-line-length-0  -g -fbounds-check -Wall -fbacktrace -ffpe-trap=invalid,zero,overflow -fsignaling-nans -fdump-core -Jobj  
FLAGS_gnu_opt   = -O3 -fdefault-real-8 -fdefault-double-8 -ffree-line-length-0 -Jobj


# intel flags
FLAGS_intel_debug = -O0 -traceback -CB -check all -warn all,nodec,interfaces -gen_interfaces
FLAGS_intel_opt   =  -r8 -O3 -ip -funroll-loops           -par-report3 -mcmodel=large -i-dynamic

# computed flags

FLAGS   = $(FLAGS_$(COMP)_$(VERSION))
LDFLAGS = $(FLAGS)
#LIBS    = -L/usr/lib64 -lfftw3



#=======================================================================
# Standard definitions
#
OFILES =  $(addprefix obj/, $(F90FILES:.f90=.o) $(F77FILES:.f=.o) $(CFILES:.c=.o))

# SGI compilers
# F90 = mpif90 -r8 -O3 -ip -funroll-loops -parallel -par-report3 -mcmodel=large -i-dynamic
#  F90 = mpif90 -r8 -O3 -ip -funroll-loops           -par-report3 -mcmodel=large -i-dynamic

#F77 = f77
#CC = cc
# LD = mpif90 -r8 -O3 -ip -funroll-loops -parallel -par-report3 -mcmodel=large -i-dynamic
#  LD = mpif90 -r8 -O3 -ip -funroll-loops           -par-report3 -mcmodel=large -i-dynamic

#VPATH = obj
#.SUFFIXES:
#.SUFFIXES: .out .o .f90 .f .c .inc .h

#=======================================================================
# Targets and dependencies
#
default:: opt
all:: msg a.out
debug:
	@make all "VERSION=debug"
prof:
	@make opt "LIBS = $(LIBS) -p"
opt:
#	@make all "FLAGS = -r8 -O3 -ip -funroll-loops -parallel -par-report3"
	@make all "VERSION=opt"
clean:
	rm -rf obj a.out *.mod
a.out: obj $(OFILES)
	@echo "$(LD) $(OFILES) $(LIBFILES) $(LIBS) "
	@$(LD) -o cfd.exe $(OFILES) $(LIBFILES) $(LIBS)
obj_cfd :  obj $(OFILES) 
	true
obj:
	[ -d obj ] || mkdir obj
msg: 
	@echo " coupler path is $(COUPLER_PATH)"

#=======================================================================
# Compilation rules
#
obj/coupler_%.o : $(COUPLER_PATH)coupler_%.f90
	$(F90) -c $(FLAGS)  $< -o $@

obj/%.o : %.f90
	$(F90) -c $(FLAGS)  $< -o $@
obj/%.o : %.f
	$(F77) $(FLAGS) -c $< -o $@
obj/%.o : %.c
	$(CC) $(FLAGS) -c $< -o $@

#dependency list
obj/fluct_random.o : obj/fluctuations_inlet_top.o
obj/main.o       :  obj/cfd_control.o 
obj/cfd_control.o  : obj/data.export_blocks.o obj/messenger.MPI.noglobal_blocks.o obj/coupler_cfd_global_data.o obj/coupler_cfd_setup.o
obj/main.o       : obj/messenger.MPI.noglobal_blocks.o
obj/coupler_cfd_setup.o : obj/data.export_blocks.o obj/simulation.o obj/coupler_cfd_global_data.o obj/mesh.export.o
obj/coupler_cfd_communication.o : obj/coupler_cfd_global_data.o 
obj/simulation.o : obj/messenger.MPI.noglobal_blocks.o obj/coupler_cfd_communication.o
obj/coupler_cfd_global_data.o : obj/data.export_blocks.o  obj/data.noglobal_blocks.o 
obj/Output_SubDomain.o : obj/messenger.MPI.noglobal_blocks.o 
obj/messenger.MPI.noglobal_blocks.o : obj/coupler_cfd_global_data.o
