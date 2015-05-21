#!/bin/bash -l

e3prep.py --job=plate_2d --do-svg

mpirun -np 5 e3mpi.exe --job=plate_2d --run

e3post.py --job=plate_2d --vtk-xml --moving-grid --tindx=all
