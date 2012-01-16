F90     := ftn
COMP    := cray
VERSION := opt

LD := $(F90)
# local libraries needed only by CFD 
ifeq ($(branch),cfd)
  ifeq ($(NAME),blayer)
    LIBFILES := $(LIBBASE)/archive.o $(LIBBASE)/library.o $(LIBBASE)/library.alc.o $(LIBBASE)/SSDfile.memory.o $(LIBBASE)/triDiagonal.o $(LIBBASE)/triDiagonal.mpi.o  $(LIBBASE)/pentaDiagonal.o $(LIBBASE)/transforms.fftw.hitachi.o $(LIBBASE)/SinCosFFT.o $(LIBBASE)/triDiagonal.mrhs.o $(LIBBASE)/triDiagonal.mpi.mrhs.o $(LIBBASE)/pentaDiagonal.mrhs.o $(LIBBASE)/pentaDiagonal.mpi.o $(LIBBASE)/polyDiagonal.o
    LIBS := -ldfftw -ldrfftw
  endif
endif


FLAGS_cray_debug_cfd :=  -e mD -J$(OBJ_DIR)
FLAGS_cray_opt_cfd   :=  -e m  -J$(OBJ_DIR)

FLAGS_cray_debug_md :=  -e mD -J$(OBJ_DIR) 
FLAGS_cray_opt_md   :=  -e m  -J$(OBJ_DIR)

# set the preprocessor flag needed to activate the coupler calls

FPP_FLAGS_no  :=  -e F -DUSE_COUPLER=0 
FPP_FLAGS_yes :=  -e F -DUSE_COUPLER=1
FPP_FLAGS     := $(FPP_FLAGS_$(USE_COUPLER))

FLAGS_$(COMP)_$(VERSION) := $(FPP_FLAGS) $(FLAGS_$(COMP)_$(VERSION)_$(branch))



