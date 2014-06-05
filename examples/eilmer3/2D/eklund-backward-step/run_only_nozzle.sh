#! /bin/bash
# run.sh
#$ -S /bin/bash
#$ -N bkstpNoz
#$ -pe orte 12
#$ -cwd
#$ -V

job=backstep
np=12

echo "Start time: "; date
e3prep.py --job=$job --zip-files
mpirun -np $np e3mpi.exe --job=$job --zip-files --run
echo "Finish time: "; date

