#! /bin/bash -I
# run.sh
#PBS -S /bin/bash
#PBS -N turbulent_plate
#PBS -l walltime=2:00:00
#PBS -l select=4:ncpus=4 -A uq-Jacobs
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
mpirun -np 16 e3mpi.exe --job=turb_flat_plate --run --max-wall-clock=7000 > LOGFILE
echo "Finish time: "; date
