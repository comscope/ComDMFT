include Makefile.in

GUTZ_ROOT2 = ${COMSUITE_BIN}


install:
	#--------------------------------------------------------
	# Build Gutzwiller solver 
	cd Gutzwiller_Solver && make && make install && cd ..
	# Copy pygtool to ${GUTZ_ROOT2}
	cp -r pygtool/* ${GUTZ_ROOT2}
	# Install pyglib
	pip install -e ./pyglib --user
	# Optionally, install local version of other python libs.
	# Install simplified version of pymatgen
	# pip install -e ./pyextern/pymatgen --user

gtool:
	cp -r pygtool/* ${GUTZ_ROOT2}

clean:
	cd Gutzwiller_Solver && make clean && cd ..
