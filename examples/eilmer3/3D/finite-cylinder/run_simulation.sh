#!/bin/bash
# run_simulation.sh
#$ -S /bin/bash
#$ -N FiniteCyl
#$ -pe orte 4
#$ -cwd
#$ -V

job=cyl
np=4

echo "Start time: "; date
mpirun -np $np e3mpi.exe --job=$job --run
# e3shared.exe --job=$job --run
echo "Finish time: "; date

