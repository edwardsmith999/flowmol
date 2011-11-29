F90     := mpif90
COMP    := intel
VERSION := opt

LD := $(F90)

# local libs, needed only for CFD 
ifeq ($(branch),cfd)
  LIBFILES := $(LIBBASE)/archive.o $(LIBBASE)/library.o $(LIBBASE)/library.alc.o $(LIBBASE)/SSDfile.memory.o $(LIBBASE)/triDiagonal.o $(LIBBASE)/triDiagonal.mpi.o $(LIBBASE)/pentaDiagonal.o $(LIBBASE)/transforms.mkl_intel.alc.o $(LIBBASE)/SinCosFFT.o $(LIBBASE)/triDiagonal.mrhs.o $(LIBBASE)/triDiagonal.mpi.mrhs.o $(LIBBASE)/pentaDiagonal.mrhs.o $(LIBBASE)/pentaDiagonal.mpi.o $(LIBBASE)/polyDiagonal.o 
endif

LIBS :=
# -warn all,nodec,interfaces
FLAGS_intel_debug_cfd := -O0 -traceback -CB -check all -module $(OBJ_DIR) 
FLAGS_intel_opt_cfd   :=  -r8 -O3 -ip -funroll-loops   -par-report3 -mcmodel=large -i-dynamic -module $(OBJ_DIR)

FLAGS_intel_debug_md := -O0 -traceback -CB -check all -module $(OBJ_DIR)
FLAGS_intel_opt_md   := -O3 -funroll-loops -module $(OBJ_DIR)

FLAGS_$(COMP)_$(VERSION) := $(FLAGS_$(COMP)_$(VERSION)_$(branch))



