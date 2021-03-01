#i&&&&&&&&&&&&&&&&&&&&&&&&&&&&        UNIX MPICH MPI       &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
include ../arch.mk

##compfl = -debug -g -CB -check bounds -traceback -check uninit -fp-model precise

wan90_dir = ../wannier90_2.1

comp = $(PF90) $(FPPFLAGS) -c $(compfl)
link = $(PF90)  $(compfl)
LIB = $(LAPACK_LIB) -L$(wan90_dir) -lwannier

#&&&&&&&&&&&&&&&&&&&& default settings &&&&&&&&&&&&&&&&&&&&

.SUFFIXES: .mod .F .o 
#compfl = -CB -check bounds -O0 -g -traceback
#compfl = -CB -check bounds -O3 -traceback
# compfl = -debug -g -CB -check bounds -traceback -check uninit -fp-model precise
#compfl = -O3
#F90 = /opt/mpich2/intel/14.0/bin/mpif90    #/opt/mpich2/intel/10.1/bin/mpif90 # ifort

#F90 = ftn
#FPPFLAGS =

#ifdef USE_HDF5
    # Turn HDF5 capabilities in the GW code on
#    FPPFLAGS += -DUSE_HDF5
    # Use the HDF5 compiler wrapper script
#    F90 = h5pfc
#F90 = /usr/common/software/hdf5-parallel/1.10.1/intel/bin/h5pfc
#endif
#
# Gfortran: The following flags are useful while developing or testing
# -DDEBUGGER: Enable code that blocks the execution until you have attached a debugger
#             You need to set
#
#                 export GW_ATTACH_DEBUGGER=Y
#
#             to block execution. Subsequently, attach debugger as 
#
#                 gdb -pid <PID>
#
#             and then in subroutine debugger
#
#                 set dosplot=.FALSE.
#
#             to resume execution.
#
# FPPFLAGS = -DMPI -DUNIX -DMDL_FIXED_FORMAT -DDEBUG -DDEBUGGER -Wuninitialized -frange-check -fbacktrace -ffpe-trap=zero,invalid,overflow
#
# The following flags are good for production versions
# FPPFLAGS = -DMPI -DUNIX
# lib1 =
# ifdef USE_HDF5
#     FPPFLAGS += -DUSE_HDF5
#     #HDF5_LIB = -L/home/hvandam/lib -lhdf5_fortran -lhdf5
# #    HDF5_LIB = -L/usr/common/software/hdf5-parallel/1.10.1/intel/lib -lhdf5_fortran -lhdf5
# #    HDF5_LIB = -L/opt/cray/pe/hdf5-parallel/1.10.0.3/INTEL/16.0/lib
# HDF5_LIB = -L/usr/common/usg/hdf5-parallel/1.10.1/intel/lib -lhdf5_fortran -lhdf5
# endif
#wan90_dir = ../wannier90_2.1
#LAPACK_LIB = 
#&&&&&&&&&&&&&&&&&&&& platform dependent settings &&&&&&&&&&&&&&&&&&&&

# GNU make does not seem to offer an if-elseif-else construct. So every
# test must appear in an if-statement of its own and it is up to you
# to ensure that there are no conflicts between the conditions.

# ifeq ($(PLATFORM),OpenSUSE-Leap-42-1-Ifort-14-0-3-OpenMPI-1-10-2)
#     F90 = /opt/ompi/intel/14.0/bin/mpif90
#     FPPFLAGS = -DMPI -DFFTMKL -DUNIX
#     LAPACK_LIB = -mkl
# endif
# #


#LIB = $(lib1) $(HDF5_LIB) $(LAPACK_LIB) -L$(wan90_dir) -lwannier
#comp = $(F90) $(FPPFLAGS) -c $(compfl)
#link = $(F90)  $(compfl)

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
o_files = \
  comwann_mod.o \
 beging.o \
   bessr.o \
   dqdall.o \
 cross.o \
 cft.o \
 diagonalization_lapack.o \
 det3.o \
 dif.o \
 fac2.o \
 fft3.o \
 fifi_j_prepare.o \
 integral_band_pw_band_mt.o \
 integral_band_pw_band.o \
 output_u_f0f2f4f6.o \
 size_shift_par.o \
 tcgstb.o \
 unport.o \
 getlimj.o \
 gbs_number.o \
 lget.o \
 prd1.o \
 prd2.o \
 radial_function.o \
 rotmat_wan_proj.o \
 rotate.o \
 sphharm.o \
 solid.o \
 sym_z_0.o \
 sym_a_0.o \
 wan_pb_rs.o \
 wannier_mmn.o \
 wannier_amn_from_win.o \
 wannier_interface.o \
 qft.o \
 init.o \
 start_par_solid.o \
 read_restart_file.o \
 set_memory.o \
 read_wan_files.o \
 read_wan_files2.o \
 read_wan_ini.o \
 

mod_files = \
 comwann_mod.mod

ifdef USE_HDF5
  o_files +=  read_wan_files_hdf5.o
  o_files +=  read_restart_file_hdf5.o
  o_files +=  write_wannier_dat_hdf5.o
  o_files +=  hdf5io_mod.o
  mod_files += hdf5io_mod.mod
endif	

TARGET = ComWann

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
	rm -rf ComWann

include ./make.depends
