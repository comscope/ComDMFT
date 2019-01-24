include ../../../Makefile.in

#FFLAGS = -O0 -check bound -check pointer -check arg_temp_created \
		-traceback


obj = GPREC.o GTIME.o GPARAM.o GCONSTANT.o GMPI.o \
	  GHDF5_BASE.o GPRIMME.o \
	  UTIL.o SPARSE.o CORRORB.o WAREHOUSE.o \
	  BANDSTRU.o DCSTD.o fermi.o \
	  GHELP.o GBROYDEN.o GKERNEL.o PSAVE.o \
	  GHYBRD.o minpack_lib.o

LIB = ../../lib/libprimme.a
INCLUDE=-I../../include

all: ${obj} CYGUTZB.o
	${PF90} -cpp ${Fmpi_define} ${FFLAGS} -o CyGutzB CYGUTZB.o \
		${obj} ${HDF5_LIB} ${MKL_LIB} ${LIB}

GPRIMME.o: ../GPRIMME.f90
	${F90} ${FIX_FORM} ${FFLAGS} -cpp -c $< ${HDF5_INCLUDE} $(INCLUDE)

%.o: ./%.f90
	${PF90} ${FREE_FORM} ${FFLAGS} -cpp ${Fmpi_define} \
			-c $< ${HDF5_INCLUDE} $(INCLUDE)

%.o: ../%.f90
	${PF90} ${FREE_FORM} ${FFLAGS} -cpp ${Fmpi_define} \
			-c $< ${HDF5_INCLUDE} $(INCLUDE)

clean:
	rm -f *o *mod *out CyGutzB

install:
	cp CyGutzB ${DESTDIR}
