F90     := ftn
COMP    := gnu
VERSION := opt

LD := $(F90)
# local libraries needed only by CFD 
ifeq ($(branch),cfd)
  ifeq ($(NAME),blayer)
    LIBFILES := $(LIBBASE)/archive.o $(LIBBASE)/library.o $(LIBBASE)/library.alc.o $(LIBBASE)/SSDfile.memory.o $(LIBBASE)/triDiagonal.o $(LIBBASE)/triDiagonal.mpi.o  $(LIBBASE)/pentaDiagonal.o $(LIBBASE)/transforms.fftw.hitachi.o $(LIBBASE)/SinCosFFT.o $(LIBBASE)/triDiagonal.mrhs.o $(LIBBASE)/triDiagonal.mpi.mrhs.o $(LIBBASE)/pentaDiagonal.mrhs.o $(LIBBASE)/pentaDiagonal.mpi.o $(LIBBASE)/polyDiagonal.o
    LIBS := -ldfftw -ldrfftw
  endif
endif


FLAGS_gnu_debug_cfd := -fdefault-real-8 -fdefault-double-8 -ffree-line-length-0  -g -fbounds-check -Wall -Wtabs -Wno-unused -fbacktrace -ffpe-trap=invalid,zero,overflow -fsignaling-nans -fdump-core -J$(OBJ_DIR)
FLAGS_gnu_opt_cfd   := -O3 -fdefault-real-8 -fdefault-double-8 -ffree-line-length-0 -J$(OBJ_DIR)

FLAGS_gnu_debug_md := -g -ffree-line-length-0 -fbounds-check -Wall -Wtabs -Wno-unused -fbacktrace -ffpe-trap=invalid,zero,overflow -fsignaling-nans -fdump-core -J$(OBJ_DIR) 
FLAGS_gnu_opt_md   := -O3 -J$(OBJ_DIR)

# set the preprocessor flag needed to activate the coupler calls

FPP_FLAGS_no  := -cpp -DUSE_COUPLER=0 
FPP_FLAGS_yes := -cpp -DUSE_COUPLER=1
FPP_FLAGS     := $(FPP_FLAGS_$(USE_COUPLER))

FLAGS_$(COMP)_$(VERSION) := $(FPP_FLAGS) $(FLAGS_$(COMP)_$(VERSION)_$(branch))



