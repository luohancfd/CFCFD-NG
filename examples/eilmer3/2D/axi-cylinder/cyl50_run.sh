#! /bin/sh
# cyl50_run.sh
e3prep.py --job=cyl50 --do-svg
time mpirun -np 4 e3mpi.exe --job=cyl50 --run

echo "At this point, we should have a new solution"
echo "Run cyl50_post.sh next"

