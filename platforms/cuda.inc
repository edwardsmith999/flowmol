CU = nvcc

FLAGS_CUDA = -O3 --ptxas-options=-v -arch=sm_11 -pg -g  

LIBS = -L /usr/local/cuda/lib64/ -lcudart
LIBSEMU	= -L /usr/local/cuda/lib64/ -lcudartemu
