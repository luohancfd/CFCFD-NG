#!/bin/bash -l
#PBS -S /bin/bash
#PBS -N optimiseM4
#PBS -q workq
#PBS -l place=scatter
#PBS -l select=1:ncpus=4:mpiprocs=4:NodeType=medium
#PBS -l walltime=168:00:00 -A uq-Jacobs
#PBS -V
#-----------------------------------------------------------

module load intel-mpi
module load intel-cc-13
module load python
echo "Where are my nodes?"
echo $PBS_NODEFILE
cat $PBS_NODEFILE
echo "-------------------------------------------"
echo "Start time: "; date
cd $PBS_O_WORKDIR
./optimise_T4_Mach_4_nozzle.py
echo "Finish time: "; date
