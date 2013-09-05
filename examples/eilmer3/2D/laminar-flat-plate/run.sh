#! /bin/bash
# run.sh
echo "Start time: "; date
mpirun -np 16 e3mpi.exe --job=lam_flat_plate --run > LOGFILE
echo "Finish time: "; date
