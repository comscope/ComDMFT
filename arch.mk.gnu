
##### fortran
F90 = mpif90  -O3  -I/usr/include/mpi -L/usr/lib/x86_64-linux-gnu -DMPI -lblas -llapack -lmpi
PF90 = mpif90 -O3  -I/usr/include/mpi -L/usr/lib/x86_64-linux-gnu -DMPI -lblas -llapack -lmpi
# compfl = -debug -g -CB -check bounds -traceback -check uninit -fp-model precise -heap-arrays
#compfl = -O3

##### f2py
fortran2python = f2py -c --fcompiler=gfortran --compiler=unix

##### phdf5
# USE_HDF5 = true
ifdef USE_HDF5
    FPPFLAGS += -DUSE_HDF5
    PF90 = h5pfc
endif

##### C and C++
CXX = cc
CXX_MPI = mpicc -DHAVE_MPI

##### lapack library
LAPACK_LIB = -llapack -lblas


##### ComCTQMC
BASE_CPPFLAGS = -DNDEBUG -I/usr/include/mpi
BASE_LIBS = -lm -llapack -lblas -lmpi
CXXFLAGS_CTQMC = -Wall -O3 -fexceptions -std=c++11 -m64

##### ComRISB
FIX_FORM = -ffixed-form
FREE_FORM = -ffree-form
#PF90 = h5pfc -I/usr/local/hdf5/include -L/usr/local/hdf5/lib
CXXFLAGS_RISB = -O2

