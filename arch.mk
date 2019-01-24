
##### fortran
F90 = ifort
PF90 = ftn
compfl = -O3

##### phdf5
USE_HDF5 = defined ### Comment out this line if you don't want to compile with hdf5 (for LDA+RISB, this line should be commented out)
ifdef USE_HDF5
    FPPFLAGS += -DUSE_HDF5
    PF90 = h5pfc
endif

##### C and C++
CXX = CC
CXX_MPI = CC

##### lapack library
LAPACK_LIB = -mkl

##### ComCTQMC
BASE_CPPFLAGS = -DNDEBUG
BASE_LIBS = -lm
CXXFLAGS_CTQMC = -std=c++11 -fexceptions -Wall -O3

##### ComRISB
FIX_FORM = -fixed
FREE_FORM = -free
PF90_RISB= h5pfc 
CXXFLAGS_RISB = -O2

