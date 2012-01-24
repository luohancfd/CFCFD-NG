#!/bin/bash -l
#PBS -S /bin/bash
#PBS -N cone20_mpi
#PBS -l select=1:ncpus=2:NodeType=medium:mpiprocs=8 -A uq-Jacobs
#PBS -l walltime=0:01:00

. /usr/share/modules/init/bash
module load python
module load intel-cc-11
module load intel-mpi
echo "Where are my nodes?"
echo $PBS_NODEFILE
cat $PBS_NODEFILE
echo "-------------------------------------------"
echo "Begin MPI job..."
date
cd $HOME/work/eilmer3/2D/cone20-simple/
mpirun -np 2 $HOME/e3bin/e3mpi.exe --job=cone20 --run > LOGFILE
echo "End MPI job."
date

