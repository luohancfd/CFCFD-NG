#!/bin/bash 

export OMP_NUM_THREADS=4

echo "part3-viscous-with-radiation"

radmodel.py -i argon-radiators-NIST-TB-QSS.py -L rad-model.lua > LOG_RADMODEL
e3prep.py --job=hemisphere > LOG_E3PREP

echo "First radiation transport calculation"
time e3rad.exe -f hemisphere -q -r > LOG_E3RAD
echo "Run e3mpi for one body length"
time mpirun -np 4 e3mpi.exe -f hemisphere -q -t 1 -r > LOG_E3MPI 

set_control_parameter.py hemisphere.control max_time 7.55855e-06 
echo "Second radiation transport calculation"
time e3rad.exe -f hemisphere -q -t 2 -r >> LOG_E3RAD
echo "Run e3mpi for one body length"
time mpirun -np 4 e3mpi.exe -f hemisphere -q -t 3 -r >> LOG_E3MPI

set_control_parameter.py hemisphere.control max_time 1.133783e-05           
echo "Third radiation transport calculation"
time e3rad.exe -f hemisphere -q -t 4 -r >> LOG_E3RAD
echo "Run e3mpi for one body length"
time mpirun -np 4 e3mpi.exe -f hemisphere -q -t 5 -r >> LOG_E3MPI

echo "Final radiation transport calculation"
time e3rad.exe -f hemisphere -q -t 6 -r >> LOG_E3RAD
