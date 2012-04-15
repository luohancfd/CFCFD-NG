#!/bin/bash 

module load openmpi-x86_64
echo "Begin MPI job..."
nohup mpirun -np 24 e3mpi.exe --job=inject --run > LOGFILE &

