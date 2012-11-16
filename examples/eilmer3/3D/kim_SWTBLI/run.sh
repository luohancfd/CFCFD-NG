job=kim_SWTBLI
np=16

e3prep.py --job=$job
mpirun -np $np e3mpi.exe --job=$job --run
