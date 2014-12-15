#! /bin/bash -I
# run.sh
#PBS -S /bin/bash
#PBS -N tc2_r64
#PBS -l walltime=72:00:00
#PBS -l select=8:ncpus=8:mpiprocs=8 -A uq-MechMinEng
#PBS -V

. /usr/share/modules/init/bash
module load intel-mpi
module load intel-cc-11
module load python

echo "Where are my nodes?"
echo $PBS_NODEFILE
cat $PBS_NODEFILE
echo "---------------------------------------------"
cd $PBS_O_WORKDIR
echo "Start time: "; date
mpirun -np 64 e3mpi.exe --job=flat_runup --run > LOGFILE
echo "Finish time: "; date
