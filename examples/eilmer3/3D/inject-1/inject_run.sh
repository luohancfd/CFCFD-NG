#! /bin/sh
# inject_run.sh
e3prep.py --job=inject
# time e3shared.exe --job=inject --run
mpirun -np 6 e3mpi.exe --job=inject --run
e3post.py --job=inject --vtk-xml
