#! /bin/bash
# run-arcus.sh
#PBS -l select=1:mpiprocs=14
#PBS -l walltime=51:00:00
#PBS -N dbl-cone
#PBS -m bea
#PBS -M peter.jacobs@eng.ox.ac.uk
#PBS -V
cd $PBS_O_WORKDIR
date
mpirun -np 14 -machinefile $PBS_NODEFILE $DATA/e3bin/e3mpi.exe \
    --job=dbl-cone --mpimap=dbl-cone.mpimap --tindx=7 --run \
    --max-wall-clock=180000 > run-arcus-3.transcript
date

echo "At this point, we should have a flow solution"
echo "Use post.sh next"

