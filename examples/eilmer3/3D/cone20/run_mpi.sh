#! /bin/sh
e3prep.py --job=cone20_3D
mpirun -np 4 e3mpi.exe --job=cone20_3D --mpimap=cone20_3D.mpimap --run
