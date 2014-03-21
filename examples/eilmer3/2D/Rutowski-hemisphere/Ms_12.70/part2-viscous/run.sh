#!/bin/bash 

echo "part2-viscous"

e3prep.py --job=hemisphere > LOG_E3PREP
echo "Adding viscous effects" 
time mpirun -np 4 e3mpi.exe -f hemisphere -q -r > LOG_E3MPI
echo "Increasing CFL number"
set_control_parameter.py hemisphere.control cfl 0.5
set_control_parameter.py hemisphere.control max_time 1.88964e-05
time mpirun -np 4 e3mpi.exe -f hemisphere -t 1 -q -r >> LOG_E3MPI

