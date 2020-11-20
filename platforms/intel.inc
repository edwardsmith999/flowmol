# Include file to be inserted in the makefile of the code to be 
# compiled in coupled mode ( or not)
#
# Provides values for the folowing make variables
# F90, LD, FPP_FLAGS, FFLAGS
#
# Needs the values of BUILD, NAME, USE_COUPLER and OBJ_DIR
#
# if BUILD is not defined it is set to opt 
#
# If NAME and USE_COUPLER are not defined warnings are printed
# if OBJ_DIR is not defined current directory is used for object files
#

# list og known builds

BUILD_LIST := debug opt profile

# check for default build

ifndef BUILD
  $(warning BUILD not define for code $(NAME), using opt)
  BUILD := opt
else
# test if the value of  BUILD is known 
  $(if $(findstring $(BUILD), $(BUILD_LIST)),,$(error unknown value of BUILD = $(BUILD)))
endif

F90   := mpif90
LD    := $(F90)

# set the preprocessor flags for coupler calls activation

ifndef USE_COUPLER
  $(warning USE_COUPLER not define for code $(NAME), setting uncoupled run)
  FPP_FLAGS     := -fpp -DUSE_COUPLER=0
else 
  FPP_FLAGS_no  := -fpp -DUSE_COUPLER=0 
  FPP_FLAGS_yes := -fpp -DUSE_COUPLER=1
  FPP_FLAGS     := $(FPP_FLAGS_$(USE_COUPLER))
endif
#
# set object files flag if OBJ_DIR is not empty
#
ifdef  OBJ_DIR
  FLAG_OBJ_DIR := -module $(OBJ_DIR)
else 
  FLAG_OBJ_DIR :=
endif
#	
# common compiler flags
#

FLAGS_BASE_debug := -g -O0 -traceback -fpp -DUSE_COUPLER=0 -fpe0 -CB $(FLAG_OBJ_DIR) #-check all -warn all -debug all -ftrapuv -fpe0
FLAGS_BASE_opt   := -O3 -funroll-loops  -fpp -DUSE_COUPLER=0 -xHost -heap-arrays 1024 -mp1 -fpconstant $(FLAG_OBJ_DIR)  #-parallel -par-report3 -mcmodel=large -i-dynamic
#FLAGS_BASE_opt   := -O3 -funroll-loops -xHost -ipo -fast -ip -heap-arrays 1024 -mp1 -fpconstant -warn general  $(FLAG_OBJ_DIR)  #-parallel -par-report3 -mcmodel=large -i-dynamic
FLAGS_BASE_profile:= -pg -g -CB -traceback -O0 -fpp -DUSE_COUPLER=0 #-O3 -funroll-loops -ip -parallel  $(FLAG_OBJ_DIR) 
#
# add  flags required by particular codes
#
ifeq ($(strip $(NAME)),coupler)

  FFLAGS := $(FLAGS_BASE_$(BUILD))  

else ifeq ($(strip $(NAME)),serial_couette)

  FFLAGS := $(FLAGS_BASE_$(BUILD)) -r8

else ifeq ($(strip $(NAME)),couette)

  LIBFILES := $(LIBBASE)/archive.o $(LIBBASE)/library.o $(LIBBASE)/library.alc.o $(LIBBASE)/SSDfile.memory.o $(LIBBASE)/triDiagonal.o $(LIBBASE)/triDiagonal.mpi.o  $(LIBBASE)/pentaDiagonal.o $(LIBBASE)/transforms.mkl_intel.alc.o $(LIBBASE)/SinCosFFT.o $(LIBBASE)/triDiagonal.mrhs.o $(LIBBASE)/triDiagonal.mpi.mrhs.o $(LIBBASE)/pentaDiagonal.mrhs.o $(LIBBASE)/pentaDiagonal.mpi.o $(LIBBASE)/polyDiagonal.o
  LIBDIR := ./../lib
  LIBFILES := $(addprefix $(LIBDIR)/,$(LIBFILES))
  #LIBS := -L/opt/intel/Compiler/11.1/046/mkl/lib/em64t/dns_lib -lmkl_em64t -lguide -lpthread
  LIBS = -mkl
  FFLAGS := -r8  $(FLAGS_BASE_$(BUILD))  

else ifeq ($(strip $(NAME)),stream_md)

  FFLAGS :=  $(FLAGS_BASE_$(BUILD))

else ifeq ($(strip $(NAME)),blayer)
#
# this must be revised
#

  LIBFILES := $(LIBBASE)/archive.o $(LIBBASE)/library.o $(LIBBASE)/library.alc.o $(LIBBASE)/SSDfile.memory.o $(LIBBASE)/triDiagonal.o $(LIBBASE)/triDiagonal.mpi.o  $(LIBBASE)/pentaDiagonal.o $(LIBBASE)/transforms.mkl_intel.alc.o $(LIBBASE)/SinCosFFT.o $(LIBBASE)/triDiagonal.mrhs.o $(LIBBASE)/triDiagonal.mpi.mrhs.o $(LIBBASE)/pentaDiagonal.mrhs.o $(LIBBASE)/pentaDiagonal.mpi.o $(LIBBASE)/polyDiagonal.o
  LIBDIR := ./../lib
  LIBFILES := $(addprefix $(LIBDIR)/,$(LIBFILES))
  #LIBS := -L/opt/intel/Compiler/11.1/046/mkl/lib/em64t/dns_lib -lmkl_em64t -lguide -lpthread
  LIBS = -mkl
  FFLAGS := $(FLAGS_BASE_$(BUILD)) -r8 -ip -parallel -par-report3 -mcmodel=large -i-dynamic

else 

    $(warning code $(NAME) not recognised, providing the base flags and preprocessor flags)
     FFLAGS := $(FLAGS_BASE_$(BUILD)) 

endif
