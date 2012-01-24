#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -N feldhuhn-mdot0
#$ -pe orte 12
echo 'Starting simulation at:'; date
mpirun -np 12 e3mpi.exe -f cone -r > outfile
echo 'Finished simulation at:'; date
