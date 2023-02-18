# ComCTQMC
ComCTQMC is a quantum impurity solver which uses the continuous time quantum Monte Carlo (CTQMC) algorithm where the action of the quantum impurity is expanded in terms of the hybridisation functions (CT-HYB). It is a stand-alone impurity solver and also embedded in dynamical mean field theory (DMFT) including ComDMFT (github.com/comscope/ComDMFT) and Portobello. 

## Features

ComCTQMC features both partition-space and worm-space sampling in order to support the measurement of all one- and two-particle Green's functions along with any static observables which can be extracted from the reduced density matrix. 

A GPU accelerated version is available for those with CUDA libraries and CUDA-capable devices (GPUs). GPUs can enable up to 600x acceleration of f-shell (14 orbital) problems or 5x acceleration of d-shell (10 orbital) problems. (Smaller problems should not use GPUs, as they will decelerate the CTQMC.)

## The Team

The project is developed by ComScope (https://www.bnl.gov/comscope/) with support from the DOE office of Basic Sciences. The main developer is Dr. Corey Melnick (cmelnick@bnl.gov).

## Citing ComCTQMC

If you use ComCTQMC, please cite our pre-print (Arxiv paper to come).

# Installation

Download the repository and unzip it, e.g., `tar -xzvf ComCTQMC.tar.gz`

The executables of ComCTQMC are compiled using `make`, and there are two executables to compile: `CTQMC` and `EVALSIM`. 

Before making these executables, one must define which compilers and libraries to use. (See the next section for a list of required libraries and compilers.) The file ComCTQMC/Makefile.in provides fields in which to set these options, along with other compiler flags you might need to get ComCTQMC working on your computer or cluster. There are a number of examples located in ComCTQMC/cluster_makefiles which can be used to compile on Cori (NERSC/LBNL), Summit (OLCF/ORNL), or a Mac. The user guide provides some additional guidance on Makefile.in

Once configured, one should execute the following commands in the ~/ComCTQMC/ directory

(cpu version) `make cpu`
(gpu version) `make gpu`

This will generate two executables: `~/ComCTQMC/bin/CTQMC` and `~/ComCTQMC/bin/EVALSIM`

If you change libraries, one should invoke `make clean` before building the executables

## Requirements

A C++11 capable compiler. The code has been tested using GNU, clang, and intel commpilers. IBM (cray) compilers are not currently supported.

BLAS and LAPACK libraries -- tested with with Intel MKL, IBM ESSL, and NETLIB-LAPACK libraries.

(optional) MPI libraries -- required for parallelization across CPUs. Tested with OpenMPI and (IBM's) Spectrum-MPI libraries.

(optional) CUDA libraries and compiler -- required for the GPU accelerated version of the code. Tested with Cuda/10.1

(optional) CUTLASS libraries -- required for the GPU accelerated version of the code. Provided with ComCTQMC.

# Usage

The general workflow required to use ComCTQMC is as follow:
1. Navigate to the working directory
2. Produce a parameter file, `params.json` (one can name this anything, provided it ends with .json)
3. Produce a file defining the hybridisation functions, `hyb.json` (one can name this anything)
4. (optional) Produce a file defining the bosonic hybridisation functions `dyn.json'
5. run the CTQMC executable
 - (mpi enabled) `mpirun -np X -npernode Y ComCTQMC/bin/CTQMC params` 
 - (otherwise) `ComCTQMC/bin/CTQMC params`
6. Run the post-processing executable
 - (mpi enabled) `mpirun -np Z -npernode Y ComCTQMC/bin/EVALSIM params`
 - (otherwise) `ComCTQMC/bin/EVALSIM params`

For a description of the input and output files, we refer the reader to the user guide UserGuide.pdf.

Examples are available in ComCTQMC/examples.
