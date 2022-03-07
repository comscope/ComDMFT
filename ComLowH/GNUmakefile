include ../arch.mk
#i&&&&&&&&&&&&&&&&&&&&&&&&&&&&        UNIX MPICH MPI       &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

#compfl = -debug -g -CB -check bounds -traceback -check uninit -fp-model precise
comp = $(PF90) $(FPPFLAGS) -c $(compfl)
link = $(PF90) $(compfl)
LIB = $(LAPACK_LIB)

TARGET = ComLowH

#&&&&&&&&&&&&&&&&&&&& default settings &&&&&&&&&&&&&&&&&&&&

.SUFFIXES: .mod .F .o 

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
o_files = \
  comlowh.o \
  comlowh_mod.o \
  pfit_mod.o \


mod_files = \
  comlowh_mod.mod \
  pfit_mod.mod \


ifdef USE_HDF5
  o_files += hdf5io_mod.o
  mod_files += hdf5io_mod.mod
  o_files +=  read_wannier_dat_hdf5.o
endif

$(TARGET): $(mod_files) $(o_files)
	$(link) $(o_files) -o $(TARGET) $(LIB)
	cp $(TARGET) ../bin
.f90.o:
	$(comp) $(INCLUDE) $<
.F.o:
	$(comp) $(INCLUDE) $<
%.mod: %.F
	$(comp) $<
%.mod: %.f90
	$(comp) $<
clean:
	rm -rf *.o
	rm -rf *.mod
	rm -rf $(TARGET)

