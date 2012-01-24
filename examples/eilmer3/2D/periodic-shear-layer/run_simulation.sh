#!/bin/bash
# run_sumulation.sh
#$ -S /bin/bash
#$ -N PeriodicShear
#$ -pe orte 4
#$ -cwd
#$ -V

job=psl
np=4

echo "Start time: "; date
mpirun -np $np e3mpi.exe --job=$job --run
# e3shared.exe --job=$job --run
echo "Finish time: "; date

