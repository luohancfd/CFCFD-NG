#! /bin/bash
# run.sh
module load openmpi-x86_64
date
mpirun -np 12 e3mpi.exe --job=cyl-flare --mpimap=cyl-flare.mpimap --run > run.transcript
date

echo "At this point, we should have a flow solution"
echo "Use post.sh next"

