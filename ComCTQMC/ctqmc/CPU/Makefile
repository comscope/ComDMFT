#include ../../include/Makefile.in
include ../../../arch.mk

BASE_LDFLAGS =

CXXFLAGS=$(CXXFLAGS_CTQMC)

CPPFLAGS = $(BASE_CPPFLAGS) $(LAPACK_CPPFLAGS)
LDFLAGS = $(BASE_LDFLAGS) $(LAPACK_LDFLAGS)
LIBS = $(BASE_LIBS) $(LAPACK_LIB)

.SUFFIXES: .C .o

HEADERS_IS = ../include/Bath.h ../include/Zahl.h ../include/Utilities.h ../include/Hyb.h ../include/Green.h ../include/Trace.h ../include/TraceAlgebra.h ../include/TraceElement.h ../include/TraceUtilities.h ../include/MarkovChain.h ../include/MonteCarlo.h ../include/Data.h ../include/Weight.h ../include/Observables.h ../../include/MPIUtilities.h ../../include/BlasLapack.h ../../include/JsonX.h ../../include/Measurements.h

o_files = \
TraceAlgebraHost.o \
ctqmc.o 

.C.o:
	$(CXX_MPI) -DHAVE_MPI $(CPPFLAGS) $(CXXFLAGS) -c $<

all:     CTQMC
	 cp CTQMC ../../../bin

CTQMC: $(o_files) $(HEADERS_IS)
	$(CXX_MPI) $(o_files) -o $@ $(LDFLAGS) $(LIBS)

clean:
	rm -f *.o CTQMC
	









 
