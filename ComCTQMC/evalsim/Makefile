include ../../arch.mk

##CUDA_CPPFLAGS = -I/hpcgpfs01/software/cuda/9.0/include/
##CUDA_LDFLAGS = -L/hpcgpfs01/software/cuda/9.0/lib64/

##NVCC = nvcc
##ARCH = -arch=sm_35
##NVCCFLAGS = -ccbin=icpc -std=c++11 -m64 -O3 
##NVCCFLAGS = -std=c++11 -m64 -O3


CPPFLAGS = $(BASE_CPPFLAGS) $(LAPACK_CPPFLAGS) $(CUDA_CPPFLAGS) 
LDFLAGS = $(BASE_LDFLAGS) $(LAPACK_LDFLAGS) $(CUDA_LDFLAGS)  
LIBS = $(BASE_LIBS) $(LAPACK_LIB) 




HEADERS = ../include/mpi/Utilities.h ../include/measurements/Measurements.h ../include/io/Vector.h ../include/io/Matrix.h ../include/JsonX.h ../include/linalg/LinAlg.h 

all:  EVALSIM 

EVALSIM:  EvalSim.C $(HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS_CTQMC) -o $@  EvalSim.C $(LDFLAGS) $(LIBS)
	cp EVALSIM ../../bin/
clean:
	rm -f *.o EVALSIM 
