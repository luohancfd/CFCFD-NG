#!/bin/bash
# run.sh
e3prep.py --job=ffs --do-svg
e3loadbalance.py --job=ffs -n 3
e3post.py --job=ffs --tindx=0 --vtk-xml

date
mpirun -np 3 e3mpi.exe --job=ffs --mpimap=ffs.mpimap --run > run.transcript
date

e3post.py --job=ffs --tindx=all --vtk-xml --add-mach

echo "At this point, we should have data to view"
