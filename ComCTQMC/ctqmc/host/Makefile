include ../../../arch.mk


##CUDA_CPPFLAGS = -I/hpcgpfs01/software/cuda/9.0/include/
##CUDA_LDFLAGS = -L/hpcgpfs01/software/cuda/9.0/lib64/

##NVCC = nvcc
##ARCH = -arch=sm_35
##NVCCFLAGS = -ccbin=icpc -std=c++11 -m64 -O3 
##NVCCFLAGS = -std=c++11 -m64 -O3


CPPFLAGS = $(BASE_CPPFLAGS) $(LAPACK_CPPFLAGS) $(CUDA_CPPFLAGS) 
LDFLAGS = $(BASE_LDFLAGS) $(LAPACK_LDFLAGS) $(CUDA_LDFLAGS)  
LIBS = $(BASE_LIBS) $(LAPACK_LIB) 


.SUFFIXES: .C .o

HEADERS_IS = ../include/bath/Bath.h ../include/Utilities.h ../include/bath/Hyb.h ../include/observables/OneParticle.h ../include/impurity/Product.h ../include/impurity/Algebra.h ../include/impurity/Node.h ../include/impurity/Operators.h ../include/MarkovChain.h ../include/MonteCarlo.h ../include/Data.h ../include/State.h ../include/Observables.h ../../include/mpi/Utilities.h ../../include/BlasLapack.h ../../include/JsonX.h ../../include/measurements/Measurements.h ../include/impurity/BitSet.h

o_files = \
AlgebraHost.o \
ctqmc.o 

.C.o:
	$(CXX_MPI) $(CPPFLAGS) $(CXXFLAGS_CTQMC) -c $<

all:     CTQMC

CTQMC: $(o_files) $(HEADERS_IS)
	$(CXX_MPI) $(o_files) -o $@ $(LDFLAGS) $(LIBS)
	cp CTQMC ../../../bin/
clean:
	rm -f *.o CTQMC
