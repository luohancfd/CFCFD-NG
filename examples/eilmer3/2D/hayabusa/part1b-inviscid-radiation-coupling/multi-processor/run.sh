#!/bin/bash

# The radiation calculation is performed with e3rad.exe, and the
# flowfield calculation with e3mpi.exe.
# This allows parallel calculation of the radiation components,
# using the previous solution (part 1).

# NOTES: 1. Set the environment variable OMP_NUM_THREADS to the 
#           the number of processors you want to use (by default
#           it should use all available processors).
#        2. Make sure e3rad has been compiled for OpenMP
#           (i.e. make e3rad TARGET=for_gnu_openmp)

python paralle_radiation_coupling.py
