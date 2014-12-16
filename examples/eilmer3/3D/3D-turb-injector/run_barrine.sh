#! /bin/bash -I
# run.sh
#PBS -S /bin/bash
#PBS -N tc3jet-96
#PBS -l walltime=150:00:00
#PBS -l select=12:ncpus=8:mpiprocs=8:mem=20g -A uq-MechMinEng
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
mpirun -np 96 e3mpi.exe --job=inject --run > LOGFILE
echo "Finish time: "; date
