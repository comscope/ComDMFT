# Quantum Mobile v21.06.04
# https://quantum-mobile.readthedocs.io/en/latest/releases/index.html

compfl = -O2 -w -fbacktrace -ffree-line-length-0

PF90 = h5pfc
F90 = h5pfc

FPPFLAGS += -DUSE_HDF5
LAPACK_LIB = -lopenblas

FIX_FORM = -ffixed-form
FREE_FORM = -ffree-form

# C/C++ compiler
CC = gcc
C++ = g++

# C compiler options.
CFLAGS = -O2
