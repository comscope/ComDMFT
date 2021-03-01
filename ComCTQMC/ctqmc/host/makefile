include ../../Makefile.in

HEADERS = ../include/bath/Bath.h ../include/Utilities.h ../include/bath/Hyb.h ../include/observables/OneParticle.h ../include/impurity/Product.h ../include/impurity/Algebra.h ../include/impurity/Node.h ../include/impurity/Operators.h ../include/MarkovChain.h ../include/MonteCarlo.h ../include/Data.h ../include/State.h ../include/Observables.h ../../include/mpi/Utilities.h ../../include/BlasLapack.h ../../include/JsonX.h ../../include/measurements/Measurements.h ../include/impurity/BitSet.h


all: CTQMC  
cpu: CTQMC

CTQMC:  ctqmc.C $(HEADERS_IS)
	$(CXX_MPI) $(CPPFLAGS) $(CXXFLAGS) -o $@  ctqmc.C $(LDFLAGS) $(LIBS)
	cp CTQMC ../../bin/.

clean:
	rm -f *.o ../../bin/CTQMC CTQMC
	









 
