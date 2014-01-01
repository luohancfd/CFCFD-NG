#!/bin/bash -l
#PBS -S /bin/bash
#PBS -N tc_flow
#PBS -q workq
#PBS -l select=5:ncpus=8:NodeType=medium:mpiprocs=8 -A uq-SCRAMSPACE
#PBS -l walltime=40:00:00
echo "-------------------------------------------"
echo "Begin MPI job..."
date
cd $PBS_O_WORKDIR
mpirun -np 40 $HOME/e3bin/e3mpi.exe --job=tc_flow_nitrogen --run --max-wall-clock=150000 > LOGFILE
echo "End MPI job."
date
