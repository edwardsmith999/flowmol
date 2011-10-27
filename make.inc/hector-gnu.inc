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


FLAGS_gnu_debug_cfd := -fdefault-real-8 -fdefault-double-8 -ffree-line-length-0  -g -fbounds-check -Wall -Wtabs -Wno-unused -fbacktrace -ffpe-trap=invalid,zero,overflow -fsignaling-nans -fdump-core -Jobj  
FLAGS_gnu_opt_cfd   := -O3 -fdefault-real-8 -fdefault-double-8 -ffree-line-length-0 -Jobj

FLAGS_gnu_debug_md := -g -fbounds-check -Wall -Wtabs -Wno-unused -fbacktrace -ffpe-trap=invalid,zero,overflow -fsignaling-nans -fdump-core -Jobj 
FLAGS_gnu_opt_md   := -O3 -Jobj

FLAGS_$(COMP)_$(VERSION) := $(FLAGS_$(COMP)_$(VERSION)_$(branch))



