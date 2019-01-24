-include ../arch.mk
-include ../../arch.mk
-include ../../../arch.mk

#=====================================================
# For Linux with intel version 11/12 on 64bit machines
#=====================================================
MPIF90=$(PF90)
COMMS=mpi

#========================================================
# Intel mkl libraries. Set LIBPATH if not in default path
#========================================================
LIBS   = $(LAPACK_LIB)
