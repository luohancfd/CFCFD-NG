#! /bin/bash -I
# run.sh
#PBS -S /bin/bash
#PBS -N turb3D
#PBS -l walltime=210:00:00
#PBS -l select=4:ncpus=4:mpiprocs=4 -A uq-MechMinEng
#PBS -V

. /usr/share/modules/init/bash
module load intel-mpi
module load intel-cc-11
module load python

echo "Where are my nodes?"
echo $PBS_NODEFILE
cat $PBS_NODEFILE
echo "----------------------------------------"
cd $PBS_O_WORKDIR
echo "Start Time: "; date
mpirun -np 16 e3mpi.exe --job=turb_flat_plate --run > LOGFILE
echo "Finish time: "; date
