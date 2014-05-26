#! /bin/bash
date
mpirun -np 4 e3mpi.exe --job=cavity --run > LOGFILE
date

