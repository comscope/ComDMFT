include ../Makefile.in

HEADERS =  

all:  EVALSIM 
cpu:  EVALSIM
gpu:  EVALSIM

EVALSIM:  Evalsim.C $(HEADERS)
	$(CXX_MPI) $(CPPFLAGS) $(CXXFLAGS) -o $@  Evalsim.C $(LDFLAGS) $(LIBS) -lmpi_cxx -lstdc++
	mv EVALSIM ../bin/.

clean:
	rm -f *.o ../bin/EVALSIM EVALSIM





 
