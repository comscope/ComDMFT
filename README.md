<!--use editor https://pandao.github.io/editor.md/en.html -->
# 1. COMSUITE 

A computational materials physics code for simulating correlated quantum materials using Dynamic Mean Field Theory (DMFT) and its extension. It can calculate the electronic structure within three different mathods:

  - charge self-consistent LDA+Gutzwiller,
  - charge self-consistent LDA+DMFT,
  - and ab initio LQSGW+DMFT
  
# 2. New version release announcement
## 2019. 1. 4
   - The first version has been released!!!
   - Please go to tutorial directory(install_directory/tutorials) to learn how to calculate the electronic structures of NiO, MnO, and FeSe. You have three choices of charge self-consistent LDA+Gutzwiller, charge self-consistent LDA+DMFT, and LQSGW+DMFT.
   
## 2020. 1. 6
   - New version released !!!
   - Now Comsuite can calculate antiferromagnetically ordered phase. Please go to tutorial directories (install_directory/tutorials/lda_dmft/NiO_afm and install_directory/tutorials/lqsgw_dmft/NiO_afm). Read pdf files to learn how to calculate the electronic structures of antiferromagnetically ordered NiO.  You have two choices of charge self-consistent LDA+DMFT and LQSGW+DMFT.

# 3. Comsuite Installation

## Prerequisites
Comsuite consists of programs, executables, and scripts, written in Fortran90, c (c++) and Python. Before you start the installation, you must make sure that the following packages are installed in your system.
  - Fortran, C, CXX compiler and blas/lapack library. The followings have been tested
    - ifort, icpc and mkl
  - MPI
    - Intel MPI (mpiifort, mpiicc, mpiicpc, mpirun, etc. Check https://software.intel.com/en-us/qualify-for-free-software)
    - open MPI (mpif90, mpicc, mpicxx, mpirun, etc. Check https://www.open-mpi.org/)
  - Python (required package : numpy, scipy, tabulate, itertools, mpi4py, cython, matplotlib, Builtins, sympy, pymatgen, pyyaml and h5py)
  
## Optional package
  - Data storage in Parallel HDF5 is also supported. Parallel HDF5 library (Check https://www.hdfgroup.org/HDF5/release/obtain5.html)

## Download COMSUITE

     git clone https://github.com/comscope/comsuite.git
    
The directory contains the following sub-directories:
-  bin -- executable binaries and scripts
- ComLowH -- program to construct low-energy Hamiltonian and tcalculate Delta
- ComWann -- program to construct Wannier function by using wannier90 package (http://wannier.org/)
- ComCoulomb -- program to calculate bosonic Weiss field
- ComCTQMC -- ctqmc impurity solver
- ComDC -- program to calculate double counted self-energy within local GW approximation
- ComRISB -- program to perform Gutzwiller-rotationally invariant slave-boson calculations.
- tutorials -- tutorials and inputfiles.
- gw -- FlapwMBPT code(https://www.bnl.gov/cmpmsd/flapwmbpt/)
- wannier90-2.1 -- the most recent version of Wannier90.

##Compile COMSUITE package.
- First, define the installation directory in the shell. For example in bash shell, use the following command adds $COMSUITE_BIN to your system $PATH

      export COMSUITE_BIN=install_directory/bin
- Then, the compilers, libraries, and flags should be defined in the arch.mk file. An example to install COMSUITE in Cori at NERSC is as follows.

       ##### fortran
       F90 = ifort
       PF90 = ftn
       compfl = -O3

       ### phdf5
       USE_HDF5 = defined  ### comment out this line if you don’t want to compile with hdf5 (for LDA+DMFT and LQSGW+DMFT)

       ifdef USE_HDF5
          FPPFLAGS += -DUSE_HDF5
	  F90 = h5pfc
       endif

       ### C and C++
       CXX = CC
       CXX_MPI = CC
	   
	    ##### lapack library
       LAPACK_LIB = -mkl

       #### ComCTQMC

       BASE_CPPFLAGS = -DNDEBUG
       BASE_LIBS = -lm

       CXXFLAGS_CTQMC = -std=c++11 -fexceptions -Wall -O3

       #### ComRISB ######################

       FIX_FORM = -fixed
       FREE_FORM = -free
       PF90_RISB=h5pfc
       CXXFLAGS_RISB = -O2
	 Below is the meaning of the each flag in the arch.mk.
  - F90 = ifort ; identify Fortran compiler for a serial fortran program
  - PF90 = ftn ; identify Fortran compiler for a MPI fortran program
  - compfl = -O3 ; compilation flag for Fortran programs
  - USE_HDF5 = defined ; compilation option to enable program to read and write data in hdf5 file format (only works for LDA+DMFT and LQSGW+DMFT for now). Comment out this line if you dont want to compile with hdf5
  - CXX = CC ; standard environment variable to identify the serial C++ compiler. 
  - CXX_MPI = CC; identify the parallel C++ compiler
  - FPPFLAGS += -DUSE_HDF5; preprocessor definition to compile with HDF5 file format
  - F90 = h5pfc; identify Fortran compiler for hdf5 file format
  - LAPACK_LIB = -mkl ; specify LAPACK library
  - BASE_CPPFLAGS = -DNDEBUG ; generic C-preprocessor flags. NDEBUG stands for no-debug code.
  - BASE_LIBS = -lm ; generic libraries that are compiler independent. The -lm flag requests the math library.
  - CXXFLAGS_CTQMC = -std=c++11 -fexceptions -Wall -O3 ; list of compilation flags for C++ compiler for ComCTQMC. The -std=c++11 specifies that the code is written based on the C++ 2011 standard. The -fexceptions flag tells the compiler that it should generate code that support exceptions. The -Wall turns on all warnings.
  - FIX_FORM = -fixed: specify Fortran standard fixed format
  - FREE_FORM = -free: specify free format
  - PF90_RISB = h5pfc: identify hdf5 supported fortran compiler for COMRISB. Here we note that ComRISB support hdf5 partially, so that “USE_HDF5=defined” should be commented out for ComRISB.
  - CXXFLAGS_RISB = -O2: c compiler option for ComRISB
- After setting up arch.mk, you need to run the following commands:

          make clean
          make

     All executable files then are in bin directory (You do not need to create bin directory by yourself). 