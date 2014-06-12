#!/bin/bash -l
#PBS -S /bin/bash
#PBS -N telbany
#PBS -q workq
#PBS -l select=12:ncpus=4:NodeType=medium:mpiprocs=4 -A uq-MechMinEng
#PBS -l walltime=80:00:00
echo "-------------------------------------------"
echo "Begin MPI job..."
date
cd $PBS_O_WORKDIR
mpirun -np 48 $HOME/e3bin/e3mpi.exe --job=couette --run -s --max-wall-clock=200000 > LOGFILE
echo "End MPI job."
date
