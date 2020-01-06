include ../../../Makefile.in

DCPP = 

FFLAGS = ${FFLAGS_RISB}

#FFLAGS = -O0 -check bound -check pointer -check arg_temp_created \
		-traceback 
#DCPP = -DDEBUG

obj = GPREC.o GTIME.o GCONSTANT.o GHDF5_BASE.o GPRIMME.o \
	  UTIL.o SPARSE.o GHDF5.o GSPCI.o GSPCI_SZ_LIB.o

LIB = ../../lib/libprimme.a
INCLUDE=-I../../include

all: ${obj} MAIN.o
	${PF90} ${FFLAGS} -o exe_spci_sz_mott MAIN.o \
		${obj} ${HDF5_LIB} ${MKL_LIB} ${LIB}

GPRIMME.o: ../GPRIMME.f90
	${F90} ${FIX_FORM} ${FFLAGS} -cpp -c $< ${HDF5_INCLUDE} $(INCLUDE)

%.o: ./%.f90
	${PF90} ${FREE_FORM} ${FFLAGS} -cpp ${DCPP} \
			-c $< ${HDF5_INCLUDE} $(INCLUDE)

%.o: ../%.f90
	${PF90} ${FREE_FORM} ${FFLAGS} -cpp  ${DCPP} \
			-c $< ${HDF5_INCLUDE} $(INCLUDE)

clean:
	rm -f *o *mod *out exe_spci_sz_mott

install:
	cp exe_spci_sz_mott ${DESTDIR}
