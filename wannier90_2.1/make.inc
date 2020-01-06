-include ../arch.mk
-include ../../arch.mk
-include ../../../arch.mk
### compfl = -debug -g -CB -check bounds -traceback -check uninit -fp-model precise
FCOPTS=$(compfl)

#=====================================================
# For Linux with intel version 11/12 on 64bit machines
#=====================================================
MPIF90=$(PF90)
COMMS=mpi

#========================================================
# Intel mkl libraries. Set LIBPATH if not in default path
#========================================================
LIBS   = $(LAPACK_LIB)
