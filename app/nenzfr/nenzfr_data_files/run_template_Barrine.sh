#!/bin/bash -l
#PBS -S /bin/bash
#PBS -N $jobName
#PBS -q workq
#PBS -l place=scatter -l select=1:ncpus=1:NodeType=fast:mpiprocs=1 -l walltime=250:00:00 -A uq-SCRAMSPACE
#PBS -V

# This file should be modified as required, however the --pe and --job options must be left as --pe=$pe
# and --job=$jobName respectively.
# Luke Doherty
# 24-02-2012

module load intel-mpi
module load intel-cc-11
module load python
echo "Where are my nodes?"
echo $PBS_NODEFILE
cat $PBS_NODEFILE
echo "-------------------------------------------"
echo "Begin Shared job..."
date

nenzfr.py --gas="air5species" --T1=300.0 --p1=250000.0 --Vs=2195.0 --pe=$pe --chem="neq" --area=1581.165 --job=$jobName --cfile="contour-t4-m10.data" --gfile="None" --exitfile="nozzle-exit.data" --nni=2700 --nnj=150 --nbi=270 --bx=1.05 --by=1.002 --max-time=0.006 --max-step=800000 > LOGFILE

echo "End Shared job."
date
