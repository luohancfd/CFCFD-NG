#!/bin/bash

# Run the simulation using e3rad.exe, not e3mpi.exe, to execute.
# This allows parallel calculation of the radiation components,
# using the previous solution (part 1).

mpirun -np 4 e3rad.exe --job=hayabusa --run
