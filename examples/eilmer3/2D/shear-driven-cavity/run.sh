#! /bin/bash
date
mpirun -np 16 e3mpi.exe --job=cavity --run > LOGFILE
date

