#!/bin/bash -l
#PBS -S /bin/bash
#PBS -N Inject16to1
#PBS -l place=scatter -l select=3:ncpus=8:NodeType=medium:mpiprocs=8 -l walltime=200:00:00 -A uq-jacobs
#PBS -V

module load python
module load intel-mpi
module load intel-cc-11
echo "Begin MPI job..."
cd $HOME/cfcfd3/examples/eilmer3/3D/rectangular-injector
mpirun -np 24 $HOME/e3bin/e3mpi.exe --job=inject --run
echo "End MPI job."
