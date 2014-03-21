#!/bin/bash 

echo "part1-inviscid"

e3prep.py --job=hemisphere > LOG_E3PREP
mpirun -np 4 e3mpi.exe -f hemisphere -r > LOG_E3MPI
