F90     := ftn
COMP    := pgi
VERISON := opt

LD := $(F90)
# local libraries needed only by CFD 
ifeq ($(branch),cfd)
  LIBFILES := $(LIBBASE)/archive.o $(LIBBASE)/library.o $(LIBBASE)/library.alc.o $(LIBBASE)/SSDfile.memory.o $(LIBBASE)/triDiagonal.o $(LIBBASE)/triDiagonal.mpi.o $(LIBBASE)/pentaDiagonal.o $(LIBBASE)/transforms.fftw.hitachi.o $(LIBBASE)/SinCosFFT.o $(LIBBASE)/triDiagonal.mrhs.o $(LIBBASE)/triDiagonal.mpi.mrhs.o $(LIBBASE)/pentaDiagonal.mrhs.o $(LIBBASE)/pentaDiagonal.mpi.o $(LIBBASE)/polyDiagonal.o
endif

LIBS := -ldfftw -ldrfftw

FLAGS_pgi_debug_cfd := -g -Mbounds -traceback -r8 -module $(OBJ_DIR)  
FLAGS_pgi_opt_cfd   := -fast -O3 -r8  -module  $(OBJ_DIR)

FLAGS_pgi_debug_md := -g -Mbounds -traceback -module $(OBJ_DIR) 
FLAGS_gnu_opt_md   := -fast -O3 -module $(OBJ_DIR)

FLAGS_$(COMP)_$(VERSION) := $(FLAGS_$(COMP)_$(VERSION)_$(branch))


