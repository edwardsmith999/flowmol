# Number of processors
NPRCS = 1

export OMP_NUM_THREADS=$(NPRCS)

# Name of the executable
PROG = SD$(NPRCS).x

# Fortran compiler
FC = ifort

# Fortran compiler options

# Debug
FFLAGS_DEBUG:= -O0 -g -traceback -CB -fpe0 -debug extended# -check all 

# Generic optimization
FFLAGS_DEFAULT:= -O3 #-axSSSE3,SSE3

# Fortran linker options
LDFLAGS =

# Directories

SRC_DIR = ./
OBJ_DIR = ./obj

# Source files
SRCS += \
$(SRC_DIR)/FEM_SD_Mod.f90 \
$(SRC_DIR)/FEM_SD_Prg.f90

# Object files
OBJS += \
$(OBJ_DIR)/FEM_SD_Mod.o \
$(OBJ_DIR)/FEM_SD_Prg.o

# Various libraries needed
MKL    =  -mkl

# MKL options for laptop, omar04, omar08
#MKL    =  -lmkl_solver -lmkl_intel -lmkl_sequential -lmkl_core -lpthread

# MKL options for cluster
#MKLPATH = ~/intel/Compiler/11.0/074/mkl/lib/em64t/
#MKL = $(MKLPATH)/libmkl_solver_lp64_sequential.a -Wl,--start-group $(MKLPATH)/libmkl_intel_lp64.a $(MKLPATH)/libmkl_sequential.a $(MKLPATH)/libmkl_core.a -Wl,--end-group -lpthread

LIBS = $(MKL)

# Rules for compilation
all: $(PROG)

$(PROG): obj setflags $(OBJS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

debug: obj changeflags $(OBJS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

setflags:
    $(eval FFLAGS = $(FFLAGS_DEFAULT))
changeflags:
    $(eval FFLAGS ?= $(FFLAGS_DEBUG))


clean:
	rm -f $(OBJS) $(OBJ_DIR)/*.mod

cleanall:
	rm -f $(PROG) $(OBJS) $(OBJ_DIR)/*.mod

obj:
	[ -d obj ] || mkdir obj

cleandat:
	rm -f ./results/*.DAT

cleanzip:
	rm -f ./results/*.gz

zip:
	gzip ./results/*.DAT

unzip:
	gunzip ./results/*.gz

# Each subdirectory must supply rules for building sources it contributes
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) -c $< -o $@ -module obj

