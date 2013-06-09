#! /bin/bash
# run-arcus.sh
#PBS -l select=1:mpiprocs=12
#PBS -l walltime=31:00:00
#PBS -N cyl-flare
#PBS -m bea
#PBS -M peter.jacobs@eng.ox.ac.uk
#PBS -V
cd $PBS_O_WORKDIR
date
mpirun -np 12 -machinefile $PBS_NODEFILE $DATA/e3bin/e3mpi.exe \
    --job=cyl-flare --mpimap=cyl-flare.mpimap --run \
    --tindx=10 --max-wall-clock=108000 > run-arcus-2.transcript
date

echo "At this point, we should have a flow solution"
echo "Use post.sh next"

