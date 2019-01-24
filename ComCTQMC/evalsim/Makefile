#include ../include/Makefile.in
include ../../arch.mk
BASE_LDFLAGS =

CXXFLAGS=$(CXXFLAGS_CTQMC)

CPPFLAGS = $(BASE_CPPFLAGS) $(LAPACK_CPPFLAGS)
LDFLAGS = $(BASE_LDFLAGS) $(LAPACK_LDFLAGS)
LIBS = $(BASE_LIBS) $(LAPACK_LIB)

HEADERS = ../include/MPIUtilities.h ../include/Measurements.h ../include/IO.h ../include/JsonX.h ../include/LinAlg.h

all:  EVALSIM
	 cp EVALSIM ../../bin
      
EVALSIM:  EvalSim.C $(HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@  EvalSim.C $(LDFLAGS) $(LIBS)

clean:
	rm -f *.o EVALSIM 





 
