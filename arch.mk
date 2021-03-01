##### fortran
F90 = ifort
PF90 = ftn
# compfl = -debug -g -CB -check bounds -traceback -check uninit -fp-model precise -heap-arrays
compfl = -O3

##### f2py
fortran2python = f2py -c --fcompiler=intelem --compiler=intelem

##### phdf5
# USE_HDF5 = true
ifdef USE_HDF5
    FPPFLAGS += -DUSE_HDF5
    PF90 = h5pfc
endif

##### C and C++
CXX = cc
CXX_MPI = CC -DHAVE_MPI

##### lapack library
LAPACK_LIB = -mkl


##### ComCTQMC
BASE_CPPFLAGS = -DNDEBUG
BASE_LIBS = -lm
CXXFLAGS_CTQMC = -Wall -O3 -fexceptions -std=c++11 -m64

##### ComRISB
FIX_FORM = -fixed
FREE_FORM = -free
PF90_risb= h5pfc
CXXFLAGS_RISB = -O2

