#!/bin/bash -l
#PBS -S /bin/bash
#PBS -N micro-combustion
#PBS -q workq
#PBS -l select=16:ncpus=4:NodeType=medium:mpiprocs=4
#PBS -A uq-MechMinEng
#PBS -l walltime=48:00:00
echo "-------------------------------------------"
echo "Begin MPI job..."
date
cd $PBS_O_WORKDIR
mpirun -np 64 e3mpi.exe --job=microchannel --run > LOGFILE
echo "End MPI job."
date
