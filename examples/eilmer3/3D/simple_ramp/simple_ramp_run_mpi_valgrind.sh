#! /bin/sh
# simple_ramp_run_mpi.sh
e3prep.py --job=simple_ramp
lamboot
mpirun -np 2 valgrind --tool=memcheck ~/e3bin/e3mpi.exe --job=simple_ramp --run --verbose
lamhalt
e3post.py --job=simple_ramp --vtk-xml --tindx=all