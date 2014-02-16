#!/bin/sh
#  back_run.sh
#  Exercise the Navier-Stokes solver for the conical nozzle 
#  as used by Back, Massier and Gier (1965) AIAA J. 3(9):1606-1614.
e3prep.py --job=back --do-svg
mpirun -np 2 e3mpi.exe --job=back --run
e3post.py --job=back --tindx=all --vtk-xml --add-mach
