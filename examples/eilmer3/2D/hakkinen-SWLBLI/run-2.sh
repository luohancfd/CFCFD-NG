#! /bin/bash
# run-2.sh
# restart the calculation where run.sh left it
module load openmpi-x86_64
date
mpirun -np 4 e3mpi.exe --job=swlbli --mpimap=swlbli.mpimap --run --tindx=10 > run-2.transcript
date

echo "At this point, we should have a flow solution"
echo "Use post.sh next"

