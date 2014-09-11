# Number of processors
NPRCS = 1

export OMP_NUM_THREADS=$(NPRCS)

# Name of the executable
PROG = SD$(NPRCS).x

# Fortran compiler
FC = ifort

# Fortran compiler options

# Debug
FFLAGS_DEBUG = -O0 -g -traceback -CB -fpe0 -debug extended -check all 

# Generic optimization
FFLAGS = -O3 #-axSSSE3,SSE3

# Fortran linker options
LDFLAGS =

# Directories

SRC_DIR = ./src
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

$(PROG): $(OBJS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

debug: $(OBJS)
	$(FC) $(FFLAGS_DEBUG) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(OBJS) $(OBJ_DIR)/*.mod

cleanall:
	rm -f $(PROG) $(OBJS) $(OBJ_DIR)/*.mod

cleandat:
	rm -f *.DAT

cleanzip:
	rm -f *.gz

zip:
	gzip *.DAT

unzip:
	gunzip *.gz

# Each subdirectory must supply rules for building sources it contributes
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) -c $< -o $@ -module obj

