#! /bin/bash -I
# run.sh
#PBS -S /bin/bash
#PBS -N tc2_fin
#PBS -l walltime=100:00:00
#PBS -l select=8:ncpus=8:mpiprocs=8:mem=20g -A uq-MechMinEng
#PBS -V

. /usr/share/modules/init/bash
module load intel-mpi
module load intel-cc-11 #was cc-13, but probably incorrect? 11 is right.
module load python

echo "Where are my nodes?"
echo $PBS_NODEFILE
cat $PBS_NODEFILE
echo "---------------------------------------------"
cd $PBS_O_WORKDIR
echo "Start time: "; date
mpirun -np 64 e3mpi.exe --job=tc2update --run > LOGFILE
#--max-wall-clock=23400
echo "Finish time: "; date
