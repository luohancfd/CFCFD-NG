#! /bin/bash
# run.sh
#PBS -S /bin/bash
#PBS -N CoaxJetsF
#PBS -l walltime=200:00:00
#PBS -l select=5:ncpus=8:mpiprocs=8 -A uq-MechMinEng
#PBS -V

job=coaxial_jets
np=40

. /usr/share/modules/init/bash
module load intel-mpi
module load intel-cc-11
module load python

echo "Where are my nodes?"
echo $PBS_NODEFILE
cat $PBS_NODEFILE
echo "----------------------------------------"
cd $PBS_O_WORKDIR
echo "Start time: "; date
mpirun -np $np e3mpi.exe --job=$job --run > LOGFILE
echo "Finish time: "; date
