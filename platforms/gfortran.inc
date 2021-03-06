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

BUILD_LIST := debug opt

# set the preprocessor flags for coupler calls activation

# check for default build
ifndef USE_COUPLER
  $(warning USE_COUPLER not define for code $(NAME), setting uncoupled run)
  FPP_FLAGS     := -cpp -DUSE_COUPLER=0
else
  FPP_FLAGS_no  := -cpp -DUSE_COUPLER=0 
  FPP_FLAGS_yes := -cpp -DUSE_COUPLER=1
  FPP_FLAGS := $(FPP_FLAGS_$(USE_COUPLER))
  $(warning USE_COUPLER defined for code $(NAME) as FPP_FLAGS_$(USE_COUPLER))
endif

# check for default build
ifndef BUILD
  $(warning BUILD not define for code $(NAME), using opt)
  BUILD := opt
else
  # test if the value of  BUILD is known 
  $(if $(findstring $(BUILD), $(BUILD_LIST)),,$(error unknown value of BUILD = $(BUILD)))
endif

ifdef NO_ASSMBLY_HEAVISIDES
  FPP_FLAGS := $(FPP_FLAGS) -DASSMBLY_HEAVISIDES=0
endif

ifdef USE_LAPACK
  FPP_FLAGS := $(FPP_FLAGS) -DUSE_LAPACK=1 -fPIC
  LDFLAGS := $(LDFLAGS) -llapack -lblas"
endif

# set the preprocessor flags for coupler calls activation

F90   := mpif90
LD    := $(F90)

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
UNDERSCORING:=-fno-underscoring

FLAGS_BASE_debug := $(FPP_FLAGS) -fbacktrace -g -ggdb -O0 -W -Wsurprising -Wampersand -Waliasing -fbounds-check -ffpe-trap=zero,overflow,underflow -J$(OBJ_DIR) $(UNDERSCORING) #-std=f2008
FLAGS_BASE_opt := -O3 -funroll-loops -J$(OBJ_DIR) $(FPP_FLAGS) $(UNDERSCORING) #-std=f2008

#
# add  flags required by particulary codes
#
ifeq ($(strip $(NAME)),coupler)

  FFLAGS := $(FLAGS_BASE_$(BUILD))

else ifeq ($(strip $(NAME)),couette)

  FFLAGS := $(FLAGS_BASE_$(BUILD)) -r8

else ifeq ($(strip $(NAME)),stream_md)

  FFLAGS :=  $(FLAGS_BASE_$(BUILD))

else ifeq ($(strip $(NAME)),blayer)
#
# this must be revised
#
  LIBFILES := $(LIBBASE)/archive.o $(LIBBASE)/library.o $(LIBBASE)/library.alc.o $(LIBBASE)/SSDfile.memory.o $(LIBBASE)/triDiagonal.o $(LIBBASE)/triDiagonal.mpi.o  $(LIBBASE)/pentaDiagonal.o $(LIBBASE)/transforms.fftw.hitachi.o $(LIBBASE)/SinCosFFT.o $(LIBBASE)/triDiagonal.mrhs.o $(LIBBASE)/triDiagonal.mpi.mrhs.o $(LIBBASE)/pentaDiagonal.mrhs.o $(LIBBASE)/pentaDiagonal.mpi.o $(LIBBASE)/polyDiagonal.o
  LIBS := -ldfftw -ldrfftw
  FFLAGS := $(FLAGS_BASE_$(BUILD)) -r8

else 

    $(warning code $(NAME) not recognised, providing the base flags and preprocessor flags)
     FFLAGS := $(FLAGS_BASE_$(BUILD)) 

endif



