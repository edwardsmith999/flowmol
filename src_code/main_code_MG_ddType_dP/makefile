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

#
# get platform to build, compiler names and flags, check the make.inc directory
#

#Check for default file
ifndef PLATFORM
    PLATFORM := $(shell cat ../../../make.inc/PLATFORM_default)
endif

ifndef PLATFORM
  platform_list := $(shell ls ../../../make.inc | grep  '.inc' | sed 's/\.inc//' )
  $(error The PLATFORM variable must be specified. Try one of the following: $(platform_list) )
else
# 
#   get the branch ( cfd or md)
#
  branch = $(if $(findstring CFD_dCSE,$(PWD)),cfd,$(if $(findstring MD_dCSE,$(PWD)),md,$(error "error in find branch")))
  include ../../../make.inc/$(PLATFORM).inc
  #Create default file for current platform
  $(shell echo $(PLATFORM) > ../../../make.inc/PLATFORM_default)
endif

# computed flags with values from PLATFORM files

FLAGS   := $(FLAGS_$(COMP)_$(VERSION))
LDFLAGS := $(FLAGS) 

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
