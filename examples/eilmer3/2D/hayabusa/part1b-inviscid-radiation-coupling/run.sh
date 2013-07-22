#!/bin/bash

# Run the simulation using e3rad.exe, not e3mpi.exe, to execute.
# This allows parallel calculation of the radiation components,
# using the previous solution (part 1).

# NOTES: 1. Set the environment variable OMP_NUM_THREADS to the 
#           the number of processors you want to use (by default
#           it should use all available processors).
#        2. Make sure e3rad has been compiled for OpenMP
#           (i.e. make e3rad TARGET=for_gnu_openmp)

e3rad.exe --job=hayabusa --run
