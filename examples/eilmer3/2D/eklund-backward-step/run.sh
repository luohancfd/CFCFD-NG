#! /bin/bash
# run.sh
#$ -S /bin/bash
#$ -N bkstp
#$ -pe orte 27
#$ -cwd
#$ -V

job=backstep
np=27

echo "Start time: "; date
e3prep.py --job=$job 
mpirun -np $np e3mpi.exe --job=$job --run 
echo "Finish time: "; date
