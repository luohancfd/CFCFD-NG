#PBS -S /bin/bash
#PBS -N hayabusa
#PBS -l place=scatter
#PBS -l select=4:ncpus=1:NodeType=medium:mpiprocs=1
#PBS -A uq-Jacobs
#PBS -l walltime=48:00:00
#PBS -V

# Seems that I have not permission to make a scratch directory...
# export SCRATCHDIR=/scratch/$USER
# if [ ! -d $SCRATCHDIR ]
# then
#     mkdir -p $SCRATCHDIR
# fi

echo "Begin mpi job..."
cd $PBS_O_WORKDIR
mpirun -np 4 e3mpi.exe -f hayabusa -r > out
echo "End mpi job."


