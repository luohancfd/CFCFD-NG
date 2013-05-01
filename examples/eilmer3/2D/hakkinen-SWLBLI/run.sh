#! /bin/bash
# run.sh
time mpirun -np 4 e3mpi.exe --job=swlbli --mpimap=swlbli.mpimap --run > run.transcript

echo "At this point, we should have a flow solution"
echo "Use post.sh next"

