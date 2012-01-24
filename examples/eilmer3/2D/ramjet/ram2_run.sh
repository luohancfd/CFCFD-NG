#! /bin/bash -I
# ram2_run.sh
#PBS -S /bin/bash
#PBS -N ram2
#PBS -l walltime=2:00:00
#PBS -l select=4:ncpus=4 -A uq-Jacobs

. /usr/share/modules/init/bash
module load intel-mpi/3.2.2.006
module load intel-cc-11
module load python

echo "Where are my nodes?"
echo $PBS_NODEFILE
cat $PBS_NODEFILE
echo "---------------------------------------------"
cd $PBS_O_WORKDIR
echo "Start time: "; date
mpirun -np 16 e3mpi.exe --job=ram2 --run --max-wall-clock=7000 > LOGFILE
echo "Finish time: "; date

