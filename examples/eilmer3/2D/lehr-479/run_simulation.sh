#!/bin/bash -l
#PBS -S /bin/bash
#PBS -N lehr
#PBS -q workq
#PBS -l select=3:ncpus=8:NodeType=medium:mpiprocs=8 -A uq-Jacobs
#PBS -l walltime=6:00:00
# Incantations to get bash to behave and the Intel MPI bits in place.
. /usr/share/modules/init/bash
module load intel-mpi/3.2.2.006
echo "Where are my nodes?"
echo $PBS_NODEFILE
cat $PBS_NODEFILE
echo "-------------------------------------------"
echo "Begin MPI job..."
date
cd $PBS_O_WORKDIR
mpirun -np 24 $HOME/e3bin/e3mpi.exe --job=lehr --run --max-wall-clock=20000 > LOGFILE
echo "End MPI job."
date
# As we leave the job, make sure that we leave no processes behind.
# (The following incantation is from Gerald Hartig.)
for i in $(cat $PBS_NODEFILE | grep -v `hostname` | sort -u); do 
    ssh $i pkill -u `whoami` 
done
killall -u `whoami` e3mpi.exe


