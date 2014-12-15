#! /bin/bash -I
# run.sh
#PBS -S /bin/bash
#PBS -N tc2_Lu_runup
#PBS -l walltime=30:00:00
#PBS -l select=4:ncpus=4:mpiprocs=4 -A uq-MechMinEng
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
mpirun -np 16 e3mpi.exe --job=flat_runup --run > LOGFILE
#--max-wall-clock=23400
echo "Finish time: "; date
