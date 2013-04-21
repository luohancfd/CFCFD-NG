# coles_run.sh
e3prep.py --job=coles --do-svg
mpirun -np 4 e3mpi.exe --job=coles --run
